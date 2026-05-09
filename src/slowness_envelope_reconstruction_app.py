# -*- coding: utf-8 -*-
"""
Lithology-constrained slowness-envelope reconstruction software.

Required input column order:
Depth, DTP, DTS, DTST, Lithology

Optional input column order:
Depth, DTP, DTS, DTST, Lithology, Profile_1, Profile_2, ...

Lithology code:
1 = Sandstone
2 = Limestone
3 = Coal
4 = Mudstone
"""

import sys
import os
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt

from matplotlib import ticker
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QFileDialog,
    QWidget,
    QLabel,
    QPushButton,
    QLineEdit,
    QTextEdit,
    QTableView,
    QProgressBar,
    QVBoxLayout,
    QHBoxLayout,
    QGridLayout,
    QTabWidget,
    QScrollArea,
    QGroupBox,
    QSplitter,
    QFrame,
    QDoubleSpinBox,
    QSpinBox,
    QCheckBox,
    QFormLayout,
    QHeaderView,
    QMessageBox,
)
from PyQt5.QtGui import QStandardItemModel, QStandardItem


plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 10


LITHOLOGY_NAMES = {
    1: "Sandstone",
    2: "Limestone",
    3: "Coal",
    4: "Mudstone",
}

LITHOLOGY_COLORS = {
    1: "yellow",
    2: "blue",
    3: "black",
    4: "green",
}

WAVE_COLORS = {
    "DTP": "#1f77b4",
    "DTS": "#d62728",
    "DTST": "#2ca02c",
}

BASELINE_COLORS = {
    "DTP": "#0b3c8c",
    "DTS": "#8b1a1a",
    "DTST": "#1b7f2a",
}

PROFILE_DEFAULT_NAMES = ["Mud", "Limestone", "Quartz", "Coal"]

PROFILE_DEFAULT_COLORS = [
    "green",
    "blue",
    "yellow",
    "black",
    "#ff7f0e",
    "#9467bd",
    "#8c564b",
    "#17becf",
    "#7f7f7f",
    "#bcbd22",
]


class NoWheelDoubleSpinBox(QDoubleSpinBox):
    def wheelEvent(self, event):
        event.ignore()


class NoWheelSpinBox(QSpinBox):
    def wheelEvent(self, event):
        event.ignore()


def read_numeric_file(file_path: str | Path) -> np.ndarray:
    file_path = Path(file_path)

    loaders = [
        lambda p: np.loadtxt(p),
        lambda p: np.loadtxt(p, delimiter="\t"),
        lambda p: np.loadtxt(p, delimiter=","),
        lambda p: np.genfromtxt(p, skip_header=1),
        lambda p: np.genfromtxt(p, delimiter="\t", skip_header=1),
        lambda p: np.genfromtxt(p, delimiter=",", skip_header=1),
    ]

    last_error = None

    for loader in loaders:
        try:
            data = loader(file_path)

            if data.ndim == 1:
                data = data.reshape(1, -1)

            if data.ndim != 2:
                continue

            data = data[~np.isnan(data).all(axis=1)]

            if data.shape[0] > 0 and data.shape[1] >= 5:
                return data

        except Exception as exc:
            last_error = exc

    raise ValueError(
        "Failed to read the input file. Required columns are: "
        "Depth, DTP, DTS, DTST, Lithology. "
        "Additional profile columns are optional."
    ) from last_error


def compute_lithology_constrained_baseline(
    depth: np.ndarray,
    slowness: np.ndarray,
    lithology: np.ndarray,
    sandstone_code: int = 1,
    window_size: int = 3,
    use_lithology_constraint: bool = True,
    clip_negative: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute the lithology-constrained floating baseline and slowness-envelope attribute.

    In sandstone intervals, local minima are used as interpolation anchors.
    In non-sandstone intervals, samples are forced to be interpolation anchors
    when lithology constraint is enabled.
    """
    depth = np.asarray(depth, dtype=float)
    slowness = np.asarray(slowness, dtype=float)
    lithology = np.asarray(lithology, dtype=int)

    if len(depth) != len(slowness) or len(depth) != len(lithology):
        raise ValueError("depth, slowness, and lithology must have the same length.")

    if len(depth) < 3:
        raise ValueError("At least three depth samples are required.")

    if window_size < 3:
        window_size = 3

    if window_size % 2 == 0:
        window_size += 1

    if window_size > len(depth):
        window_size = len(depth) if len(depth) % 2 == 1 else len(depth) - 1

    half_window = window_size // 2
    anchor_mask = np.zeros(len(depth), dtype=bool)

    anchor_mask[0] = True
    anchor_mask[-1] = True

    for i in range(half_window, len(depth) - half_window):
        if use_lithology_constraint and lithology[i] != sandstone_code:
            anchor_mask[i] = True
            continue

        window = slowness[i - half_window:i + half_window + 1]
        center_value = slowness[i]

        if center_value == np.min(window) and np.sum(window == center_value) == 1:
            anchor_mask[i] = True

    if use_lithology_constraint:
        anchor_mask[lithology != sandstone_code] = True

    anchor_indices = np.where(anchor_mask)[0]

    if len(anchor_indices) < 2:
        anchor_indices = np.array([0, len(depth) - 1], dtype=int)

    baseline = np.interp(
        depth,
        depth[anchor_indices],
        slowness[anchor_indices],
    )

    envelope = slowness - baseline

    if clip_negative:
        envelope = np.maximum(envelope, 0.0)

    return baseline, envelope, anchor_indices


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Lithology-constrained Slowness-envelope Reconstruction")
        self.resize(1520, 940)

        self.data = None
        self.depth = None

        self.dtp = None
        self.dts = None
        self.dtst = None

        self.lithology = None
        self.profiles = None
        self.profile_names = []

        self.dtp_baseline = None
        self.dts_baseline = None
        self.dtst_baseline = None

        self.s_dtp = None
        self.s_dts = None
        self.s_dtst = None

        self.anchor_dtp = None
        self.anchor_dts = None
        self.anchor_dtst = None

        self.segment_id_result = None

        self.log_counter = 0
        self.segment_cards = []

        self.canvas = None
        self.toolbar = None
        self.current_figure = None
        self.depth_scrollbar = None

        self.full_depth_min = None
        self.full_depth_max = None
        self.current_view_min = None
        self.current_view_max = None

        self._build_ui()
        self._apply_styles()
        self._set_global_font()

    def _set_global_font(self):
        app_font = QtGui.QFont("Arial")
        app_font.setPointSize(10)
        QApplication.instance().setFont(app_font)

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)

        main_layout = QVBoxLayout(central)
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(8)

        top_card = QFrame()
        top_card.setObjectName("TopCard")
        top_layout = QHBoxLayout(top_card)
        top_layout.setContentsMargins(12, 12, 12, 12)
        top_layout.setSpacing(10)

        self.file_edit = QLineEdit()
        self.file_edit.setPlaceholderText("Select input TXT, CSV, or DAT file")
        self.file_edit.setMinimumHeight(38)

        self.btn_select = QPushButton("Select File")
        self.btn_read = QPushButton("Read Data")
        self.btn_run = QPushButton("Run")
        self.btn_export = QPushButton("Export Results")
        self.btn_save_figure = QPushButton("Export Figure")

        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setFixedWidth(220)

        top_layout.addWidget(self.file_edit, 1)
        top_layout.addWidget(self.btn_select)
        top_layout.addWidget(self.btn_read)
        top_layout.addWidget(self.btn_run)
        top_layout.addWidget(self.btn_export)
        top_layout.addWidget(self.btn_save_figure)
        top_layout.addWidget(self.progress_bar)

        main_layout.addWidget(top_card)

        self.tabs = QTabWidget()
        self.tabs.setTabPosition(QTabWidget.North)
        main_layout.addWidget(self.tabs, 1)

        self._build_page_data()
        self._build_page_analysis()

        self.tabs.addTab(self.page_data, "Data Import and Preview")
        self.tabs.addTab(self.page_analysis, "Parameter Settings and Plotting")

        footer = QFrame()
        footer.setObjectName("FooterCard")
        footer_layout = QHBoxLayout(footer)
        footer_layout.setContentsMargins(10, 6, 10, 6)

        self.footer_left = QLineEdit("Module: lithology-constrained slowness-envelope reconstruction")
        self.footer_mid = QLineEdit("Array acoustic logging")
        self.footer_right = QLineEdit("DTP / DTS / DTST")

        for widget in [self.footer_left, self.footer_mid, self.footer_right]:
            widget.setReadOnly(True)
            widget.setFrame(False)

        self.footer_mid.setAlignment(QtCore.Qt.AlignCenter)
        self.footer_right.setAlignment(QtCore.Qt.AlignRight)

        footer_layout.addWidget(self.footer_left, 1)
        footer_layout.addWidget(self.footer_mid, 1)
        footer_layout.addWidget(self.footer_right, 1)

        main_layout.addWidget(footer)

        self.btn_select.clicked.connect(self.on_select_file)
        self.btn_read.clicked.connect(self.on_read_data)
        self.btn_run.clicked.connect(self.on_run)
        self.btn_export.clicked.connect(self.out)
        self.btn_save_figure.clicked.connect(self.save_figure)
        self.btn_add_segment.clicked.connect(self.add_segment)
        self.btn_remove_segment.clicked.connect(self.remove_segment)

    def _build_page_data(self):
        self.page_data = QWidget()
        layout = QVBoxLayout(self.page_data)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(8)

        info_frame = QFrame()
        info_frame.setObjectName("InfoCard")
        info_layout = QVBoxLayout(info_frame)
        info_layout.setContentsMargins(12, 10, 12, 10)

        self.data_info_label = QLabel(
            "No data loaded.\n"
            "Input column order: Depth, DTP, DTS, DTST, Lithology, optional profile columns.\n"
            "Lithology code: 1 = Sandstone, 2 = Limestone, 3 = Coal, 4 = Mudstone."
        )
        self.data_info_label.setWordWrap(True)
        info_layout.addWidget(self.data_info_label)

        layout.addWidget(info_frame)

        splitter = QSplitter(QtCore.Qt.Horizontal)

        self.table_preview = QTableView()
        self.table_preview.setAlternatingRowColors(True)
        self.table_preview.setSortingEnabled(False)
        self.table_preview.verticalHeader().setVisible(False)
        self.table_preview.horizontalHeader().setStretchLastSection(True)
        self.table_preview.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        splitter.addWidget(self.table_preview)

        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        splitter.addWidget(self.log_text)

        splitter.setSizes([980, 420])
        layout.addWidget(splitter, 1)

    def _build_page_analysis(self):
        self.page_analysis = QWidget()
        page_layout = QHBoxLayout(self.page_analysis)
        page_layout.setContentsMargins(8, 8, 8, 8)
        page_layout.setSpacing(8)

        splitter = QSplitter(QtCore.Qt.Horizontal)
        page_layout.addWidget(splitter)

        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        left_layout.setContentsMargins(6, 6, 6, 6)
        left_layout.setSpacing(8)

        segment_ctrl_frame = QFrame()
        segment_ctrl_frame.setObjectName("SegmentControlCard")
        segment_ctrl_layout = QHBoxLayout(segment_ctrl_frame)
        segment_ctrl_layout.setContentsMargins(10, 10, 10, 10)

        segment_label = QLabel("Segment Settings")
        segment_label.setObjectName("SectionTitle")

        self.btn_add_segment = QPushButton("+ Add Segment")
        self.btn_remove_segment = QPushButton("- Remove Last Segment")

        segment_ctrl_layout.addWidget(segment_label)
        segment_ctrl_layout.addStretch()
        segment_ctrl_layout.addWidget(self.btn_add_segment)
        segment_ctrl_layout.addWidget(self.btn_remove_segment)

        left_layout.addWidget(segment_ctrl_frame)

        self.segment_scroll = QScrollArea()
        self.segment_scroll.setWidgetResizable(True)
        self.segment_scroll.setFrameShape(QtWidgets.QFrame.NoFrame)

        self.segment_container = QWidget()
        self.segment_container_layout = QVBoxLayout(self.segment_container)
        self.segment_container_layout.setContentsMargins(4, 4, 4, 4)
        self.segment_container_layout.setSpacing(10)
        self.segment_container_layout.addStretch()

        self.segment_scroll.setWidget(self.segment_container)
        left_layout.addWidget(self.segment_scroll, 1)

        log_frame2 = QGroupBox("Run Log")
        log_layout2 = QVBoxLayout(log_frame2)
        log_layout2.setContentsMargins(8, 8, 8, 8)

        self.log_text_page2 = QTextEdit()
        self.log_text_page2.setReadOnly(True)
        self.log_text_page2.setMinimumHeight(210)
        log_layout2.addWidget(self.log_text_page2)

        left_layout.addWidget(log_frame2)

        splitter.addWidget(left_panel)

        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        right_layout.setContentsMargins(6, 6, 6, 6)
        right_layout.setSpacing(8)

        self.plot_frame = QFrame()
        self.plot_frame.setObjectName("PlotCard")
        self.plot_frame.setLayout(QVBoxLayout())
        self.plot_frame.layout().setContentsMargins(0, 0, 0, 0)
        right_layout.addWidget(self.plot_frame, 1)

        splitter.addWidget(right_panel)
        splitter.setSizes([480, 1040])

    def _apply_styles(self):
        self.setStyleSheet("""
        * {
            font-family: "Arial";
            font-size: 10pt;
        }

        QMainWindow {
            background: #f2f4f7;
        }

        QFrame#TopCard, QFrame#InfoCard, QFrame#FooterCard, QFrame#SegmentControlCard, QFrame#PlotCard {
            background: #ffffff;
            border: 1px solid #d9dde3;
            border-radius: 10px;
        }

        QFrame {
            background: #ffffff;
            border: 1px solid #e0e4ea;
            border-radius: 8px;
        }

        QLabel#SectionTitle {
            font-size: 11pt;
            font-weight: bold;
            color: #2f3a45;
        }

        QLineEdit, QTextEdit, QTableView, QDoubleSpinBox, QSpinBox {
            background: #ffffff;
            border: 1px solid #cfd5dc;
            border-radius: 6px;
            padding: 4px 6px;
            selection-background-color: #d7dde5;
        }

        QPushButton {
            background: #e5e7eb;
            border: 1px solid #bcc3cb;
            border-radius: 7px;
            padding: 8px 14px;
            min-height: 18px;
            color: #111827;
            font-weight: 600;
        }

        QPushButton:hover {
            background: #d9dde2;
            border: 1px solid #aab2bc;
            color: #0f172a;
        }

        QPushButton:pressed {
            background: #cfd4da;
            color: #0b1220;
        }

        QProgressBar {
            background: #eceff3;
            border: 1px solid #cfd5dc;
            border-radius: 6px;
            text-align: center;
            color: #2f3a45;
        }

        QProgressBar::chunk {
            background: #22c55e;
            border-radius: 5px;
        }

        QTabWidget::pane {
            border: 1px solid #d9dde3;
            background: #ffffff;
            border-radius: 8px;
        }

        QTabBar::tab {
            background: #e6e9ed;
            color: #3a4652;
            padding: 8px 18px;
            margin-right: 2px;
            border-top-left-radius: 8px;
            border-top-right-radius: 8px;
        }

        QTabBar::tab:selected {
            background: #ffffff;
            font-weight: bold;
            color: #2f3a45;
        }

        QGroupBox {
            font-weight: bold;
            border: 1px solid #d9dde3;
            border-radius: 8px;
            margin-top: 10px;
            background: #fafbfc;
        }

        QGroupBox::title {
            subcontrol-origin: margin;
            left: 10px;
            padding: 0 4px;
            color: #2f3a45;
        }

        QHeaderView::section {
            background: #eceff3;
            color: #2f3a45;
            border: 1px solid #d7dce3;
            padding: 5px;
            font-weight: bold;
        }

        QScrollBar:vertical {
            background: #eef1f4;
            width: 14px;
            border-radius: 6px;
            margin: 0;
        }

        QScrollBar::handle:vertical {
            background: #9aa4af;
            min-height: 40px;
            border-radius: 6px;
        }

        QScrollBar::handle:vertical:hover {
            background: #7f8b97;
        }

        QCheckBox {
            spacing: 6px;
        }
        """)

    def log_message(self, message: str, level: str = "INFO"):
        self.log_counter += 1
        timestamp = QtCore.QDateTime.currentDateTime().toString("yyyy-MM-dd HH:mm:ss")

        level = level.upper().strip()
        color = {
            "ERROR": "#b42318",
            "WARNING": "#b54708",
            "SUCCESS": "#027a48",
            "INFO": "#344054",
        }.get(level, "#344054")

        html = (
            f"<div style='margin:4px 0;'>"
            f"<span style='color:#667085;'>[{self.log_counter:03d}] {timestamp}</span> "
            f"<span style='color:{color}; font-weight:600;'>[{level}]</span> "
            f"<span style='color:#1f2937;'>{message}</span>"
            f"</div>"
        )

        self.log_text.append(html)
        self.log_text.ensureCursorVisible()
        self.log_text_page2.append(html)
        self.log_text_page2.ensureCursorVisible()

    def on_select_file(self):
        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select File",
            "",
            "Data Files (*.txt *.csv *.dat);;Text Files (*.txt);;CSV Files (*.csv);;All Files (*)",
            options=options,
        )
        if file_path:
            self.file_edit.setText(file_path)
            self.log_message(f"Input file selected: {file_path}", "INFO")

    def on_read_data(self):
        filename = self.file_edit.text().strip()

        if not filename:
            self.log_message("No input file selected.", "WARNING")
            return

        if not os.path.exists(filename):
            self.log_message("The selected file does not exist.", "ERROR")
            return

        try:
            self.progress_bar.setValue(10)
            self.log_message("Reading input data file.", "INFO")

            data = read_numeric_file(filename)

            if data.shape[1] < 5:
                raise ValueError(
                    "The input file must contain at least 5 columns: "
                    "Depth, DTP, DTS, DTST, Lithology."
                )

            self.data = data
            self.depth = data[:, 0]
            self.dtp = data[:, 1]
            self.dts = data[:, 2]
            self.dtst = data[:, 3]
            self.lithology = np.rint(data[:, 4]).astype(int)

            if data.shape[1] > 5:
                self.profiles = data[:, 5:]
            else:
                self.profiles = np.empty((len(self.depth), 0), dtype=float)

            self.profile_names = []
            for index in range(self.profiles.shape[1]):
                if index < len(PROFILE_DEFAULT_NAMES):
                    self.profile_names.append(PROFILE_DEFAULT_NAMES[index])
                else:
                    self.profile_names.append(f"Profile_{index + 1}")

            invalid_codes = set(np.unique(self.lithology).tolist()) - set(LITHOLOGY_NAMES.keys())
            if invalid_codes:
                raise ValueError(
                    f"Invalid lithology code detected: {sorted(invalid_codes)}. "
                    "Allowed values are 1, 2, 3, and 4."
                )

            if np.any(np.diff(self.depth) == 0):
                raise ValueError("Duplicate depth values were detected.")

            self.progress_bar.setValue(45)
            self.show_data_in_table()

            self.data_info_label.setText(
                f"Data loaded: rows = {self.data.shape[0]}, columns = {self.data.shape[1]}; "
                f"depth range = {self.depth.min():.3f} to {self.depth.max():.3f}\n"
                "Input column order: Depth, DTP, DTS, DTST, Lithology, optional profile columns.\n"
                "Lithology code: 1 = Sandstone, 2 = Limestone, 3 = Coal, 4 = Mudstone."
            )

            self.initialize_default_segments()

            self.progress_bar.setValue(100)
            self.log_message("Data loading completed.", "SUCCESS")

            if self.profiles.shape[1] > 0:
                self.log_message(f"Profile columns detected: {', '.join(self.profile_names)}.", "INFO")
            else:
                self.log_message("No profile columns were detected.", "INFO")

        except Exception as exc:
            self.progress_bar.setValue(0)
            self.log_message(f"Failed to read data: {exc}", "ERROR")
            QMessageBox.critical(self, "Read Error", str(exc))

    def show_data_in_table(self):
        model = QStandardItemModel()

        headers = [
            "Depth (m)",
            "DTP (us/ft)",
            "DTS (us/ft)",
            "DTST (us/ft)",
            "Lithology",
        ]

        if self.profiles is not None and self.profiles.shape[1] > 0:
            headers += self.profile_names

        model.setHorizontalHeaderLabels(headers)

        preview_data = np.column_stack([
            self.depth,
            self.dtp,
            self.dts,
            self.dtst,
            self.lithology,
        ])

        if self.profiles is not None and self.profiles.shape[1] > 0:
            preview_data = np.column_stack([preview_data, self.profiles])

        max_rows = min(500, preview_data.shape[0])

        for row in preview_data[:max_rows]:
            items = [QStandardItem(str(value)) for value in row]
            model.appendRow(items)

        self.table_preview.setModel(model)
        self.table_preview.resizeColumnsToContents()

    def initialize_default_segments(self):
        for card in self.segment_cards:
            card["group"].setParent(None)

        self.segment_cards.clear()

        if self.depth is None or len(self.depth) == 0:
            return

        self.add_segment(
            top=float(np.min(self.depth)),
            bottom=float(np.max(self.depth)),
            copy_from_last=False,
        )

    def add_segment(self, checked=False, top=None, bottom=None, copy_from_last=True):
        if self.depth is None:
            self.log_message("Load data before adding a segment.", "WARNING")
            return

        seg_index = len(self.segment_cards) + 1

        if top is None:
            top = float(np.min(self.depth))

        if bottom is None:
            bottom = float(np.max(self.depth))

        group = QGroupBox(f"Segment {seg_index}")
        group_layout = QVBoxLayout(group)
        group_layout.setContentsMargins(10, 12, 10, 10)
        group_layout.setSpacing(8)

        depth_box = QGroupBox("Depth Range")
        depth_form = QFormLayout(depth_box)

        top_spin = self._make_double_spinbox(top, 3)
        bottom_spin = self._make_double_spinbox(bottom, 3)

        depth_form.addRow("Top depth", top_spin)
        depth_form.addRow("Bottom depth", bottom_spin)
        group_layout.addWidget(depth_box)

        algorithm_box = QGroupBox("Algorithm Parameters")
        algorithm_form = QFormLayout(algorithm_box)

        widgets = {}

        widgets["window_size"] = NoWheelSpinBox()
        widgets["window_size"].setRange(3, 99)
        widgets["window_size"].setSingleStep(2)
        widgets["window_size"].setValue(3)

        widgets["use_lithology_constraint"] = QCheckBox("Enable lithology constraint")
        widgets["use_lithology_constraint"].setChecked(True)

        widgets["clip_negative"] = QCheckBox("Clip negative values to zero")
        widgets["clip_negative"].setChecked(True)

        widgets["show_anchors"] = QCheckBox("Show interpolation anchors")
        widgets["show_anchors"].setChecked(False)

        algorithm_form.addRow("Sliding window length", widgets["window_size"])
        algorithm_form.addRow(widgets["use_lithology_constraint"])
        algorithm_form.addRow(widgets["clip_negative"])
        algorithm_form.addRow(widgets["show_anchors"])

        group_layout.addWidget(algorithm_box)

        display_box = QGroupBox("Display Parameters")
        display_grid = QGridLayout(display_box)

        widgets["s_dtp_threshold"] = self._make_double_spinbox(0.0, 3)
        widgets["s_dts_threshold"] = self._make_double_spinbox(0.0, 3)
        widgets["s_dtst_threshold"] = self._make_double_spinbox(0.0, 3)

        widgets["dtp_min"] = self._make_double_spinbox(40.0, 2)
        widgets["dtp_max"] = self._make_double_spinbox(140.0, 2)
        widgets["dts_min"] = self._make_double_spinbox(80.0, 2)
        widgets["dts_max"] = self._make_double_spinbox(300.0, 2)
        widgets["dtst_min"] = self._make_double_spinbox(200.0, 2)
        widgets["dtst_max"] = self._make_double_spinbox(300.0, 2)

        display_grid.addWidget(QLabel("S_DTP threshold"), 0, 0)
        display_grid.addWidget(widgets["s_dtp_threshold"], 0, 1)
        display_grid.addWidget(QLabel("S_DTS threshold"), 1, 0)
        display_grid.addWidget(widgets["s_dts_threshold"], 1, 1)
        display_grid.addWidget(QLabel("S_DTST threshold"), 2, 0)
        display_grid.addWidget(widgets["s_dtst_threshold"], 2, 1)

        display_grid.addWidget(QLabel("DTP min"), 3, 0)
        display_grid.addWidget(widgets["dtp_min"], 3, 1)
        display_grid.addWidget(QLabel("DTP max"), 4, 0)
        display_grid.addWidget(widgets["dtp_max"], 4, 1)
        display_grid.addWidget(QLabel("DTS min"), 5, 0)
        display_grid.addWidget(widgets["dts_min"], 5, 1)
        display_grid.addWidget(QLabel("DTS max"), 6, 0)
        display_grid.addWidget(widgets["dts_max"], 6, 1)
        display_grid.addWidget(QLabel("DTST min"), 7, 0)
        display_grid.addWidget(widgets["dtst_min"], 7, 1)
        display_grid.addWidget(QLabel("DTST max"), 8, 0)
        display_grid.addWidget(widgets["dtst_max"], 8, 1)

        group_layout.addWidget(display_box)

        if copy_from_last and self.segment_cards:
            self._copy_segment_parameters(self.segment_cards[-1]["widgets"], widgets)

        self.segment_container_layout.insertWidget(self.segment_container_layout.count() - 1, group)

        self.segment_cards.append({
            "group": group,
            "top_spin": top_spin,
            "bottom_spin": bottom_spin,
            "widgets": widgets,
        })

        self.log_message(f"Segment {seg_index} added.", "INFO")

    def remove_segment(self):
        if not self.segment_cards:
            return

        if len(self.segment_cards) == 1:
            self.log_message("At least one segment must be retained.", "WARNING")
            return

        card = self.segment_cards.pop()
        card["group"].setParent(None)
        self._refresh_segment_titles()
        self.log_message("The last segment was removed.", "INFO")

    def _refresh_segment_titles(self):
        for index, card in enumerate(self.segment_cards, start=1):
            card["group"].setTitle(f"Segment {index}")

    def _copy_segment_parameters(self, src, dst):
        for key in dst:
            if key not in src:
                continue

            if isinstance(dst[key], QCheckBox):
                dst[key].setChecked(src[key].isChecked())
            elif isinstance(dst[key], (QDoubleSpinBox, QSpinBox)):
                dst[key].setValue(src[key].value())

    def _make_double_spinbox(self, value=0.0, decimals=2):
        spin_box = NoWheelDoubleSpinBox()
        spin_box.setDecimals(decimals)
        spin_box.setRange(-1e12, 1e12)
        spin_box.setSingleStep(0.1)
        spin_box.setValue(value)
        spin_box.setMinimumHeight(28)
        return spin_box

    def collect_segment_configs(self):
        configs = []

        for index, card in enumerate(self.segment_cards, start=1):
            widgets = card["widgets"]

            window_size = int(widgets["window_size"].value())
            if window_size % 2 == 0:
                window_size += 1
                widgets["window_size"].setValue(window_size)

            configs.append({
                "segment_id": index,
                "top": float(card["top_spin"].value()),
                "bottom": float(card["bottom_spin"].value()),
                "window_size": window_size,
                "use_lithology_constraint": widgets["use_lithology_constraint"].isChecked(),
                "clip_negative": widgets["clip_negative"].isChecked(),
                "show_anchors": widgets["show_anchors"].isChecked(),
                "s_dtp_threshold": float(widgets["s_dtp_threshold"].value()),
                "s_dts_threshold": float(widgets["s_dts_threshold"].value()),
                "s_dtst_threshold": float(widgets["s_dtst_threshold"].value()),
                "dtp_min": float(widgets["dtp_min"].value()),
                "dtp_max": float(widgets["dtp_max"].value()),
                "dts_min": float(widgets["dts_min"].value()),
                "dts_max": float(widgets["dts_max"].value()),
                "dtst_min": float(widgets["dtst_min"].value()),
                "dtst_max": float(widgets["dtst_max"].value()),
            })

        return configs

    def on_run(self):
        if self.data is None:
            self.log_message("No data are available. Please read data first.", "WARNING")
            return

        try:
            self.progress_bar.setValue(10)

            configs = self.collect_segment_configs()
            self._validate_segments(configs)

            self.progress_bar.setValue(25)
            self.log_message("Running lithology-constrained slowness-curve reconstruction.", "INFO")
            self.calculate_results(configs)

            self.progress_bar.setValue(75)
            self.log_message("Generating the integrated logging plot.", "INFO")
            self.plot_graph(configs)

            self.progress_bar.setValue(100)
            self.tabs.setCurrentIndex(1)
            self.log_message("Processing completed.", "SUCCESS")

        except Exception as exc:
            self.progress_bar.setValue(0)
            self.log_message(f"Processing failed: {exc}", "ERROR")
            QMessageBox.critical(self, "Run Error", str(exc))

    def _validate_segments(self, configs):
        depth_min = float(np.min(self.depth))
        depth_max = float(np.max(self.depth))

        for cfg in configs:
            top = min(cfg["top"], cfg["bottom"])
            bottom = max(cfg["top"], cfg["bottom"])

            if bottom < depth_min or top > depth_max:
                raise ValueError(f"Segment {cfg['segment_id']} is outside the data depth range.")

        sorted_configs = sorted(configs, key=lambda item: min(item["top"], item["bottom"]))

        for index in range(len(sorted_configs) - 1):
            first_bottom = max(sorted_configs[index]["top"], sorted_configs[index]["bottom"])
            second_top = min(sorted_configs[index + 1]["top"], sorted_configs[index + 1]["bottom"])

            if second_top < first_bottom:
                self.log_message(
                    f"Segment {sorted_configs[index]['segment_id']} overlaps with "
                    f"Segment {sorted_configs[index + 1]['segment_id']}. "
                    "Results in the overlapping interval will be overwritten by the later segment.",
                    "WARNING",
                )

    def calculate_results(self, configs):
        sample_count = len(self.depth)

        self.dtp_baseline = np.full(sample_count, np.nan)
        self.dts_baseline = np.full(sample_count, np.nan)
        self.dtst_baseline = np.full(sample_count, np.nan)

        self.s_dtp = np.full(sample_count, np.nan)
        self.s_dts = np.full(sample_count, np.nan)
        self.s_dtst = np.full(sample_count, np.nan)

        self.anchor_dtp = np.zeros(sample_count, dtype=bool)
        self.anchor_dts = np.zeros(sample_count, dtype=bool)
        self.anchor_dtst = np.zeros(sample_count, dtype=bool)

        self.segment_id_result = np.full(sample_count, -1, dtype=int)

        for cfg in configs:
            self._calculate_single_segment(cfg)

        self.log_message("All segments were processed and merged successfully.", "SUCCESS")

    def _calculate_single_segment(self, cfg):
        segment_id = cfg["segment_id"]
        top = min(cfg["top"], cfg["bottom"])
        bottom = max(cfg["top"], cfg["bottom"])

        mask = (self.depth >= top) & (self.depth <= bottom)
        indices = np.where(mask)[0]

        if len(indices) < 3:
            self.log_message(f"Segment {segment_id} has too few samples and was skipped.", "WARNING")
            return

        depth_segment = self.depth[indices]
        lithology_segment = self.lithology[indices]

        use_lithology_constraint = cfg["use_lithology_constraint"]
        clip_negative = cfg["clip_negative"]
        window_size = cfg["window_size"]

        dtp_baseline, s_dtp, anchor_dtp_index = compute_lithology_constrained_baseline(
            depth_segment,
            self.dtp[indices],
            lithology_segment,
            sandstone_code=1,
            window_size=window_size,
            use_lithology_constraint=use_lithology_constraint,
            clip_negative=clip_negative,
        )

        dts_baseline, s_dts, anchor_dts_index = compute_lithology_constrained_baseline(
            depth_segment,
            self.dts[indices],
            lithology_segment,
            sandstone_code=1,
            window_size=window_size,
            use_lithology_constraint=use_lithology_constraint,
            clip_negative=clip_negative,
        )

        dtst_baseline, s_dtst, anchor_dtst_index = compute_lithology_constrained_baseline(
            depth_segment,
            self.dtst[indices],
            lithology_segment,
            sandstone_code=1,
            window_size=window_size,
            use_lithology_constraint=use_lithology_constraint,
            clip_negative=clip_negative,
        )

        self.dtp_baseline[indices] = dtp_baseline
        self.dts_baseline[indices] = dts_baseline
        self.dtst_baseline[indices] = dtst_baseline

        self.s_dtp[indices] = s_dtp
        self.s_dts[indices] = s_dts
        self.s_dtst[indices] = s_dtst

        self.anchor_dtp[indices[anchor_dtp_index]] = True
        self.anchor_dts[indices[anchor_dts_index]] = True
        self.anchor_dtst[indices[anchor_dtst_index]] = True

        self.segment_id_result[indices] = segment_id

        self.log_message(
            f"Segment {segment_id} processed: depth {top:.3f} to {bottom:.3f}, "
            f"samples = {len(indices)}.",
            "INFO",
        )

    def plot_graph(self, configs):
        if self.depth is None:
            return

        fig, axes = plt.subplots(1, 8, figsize=(18.8, 8.4), sharey=True)
        plt.subplots_adjust(left=0.05, right=0.985, top=0.94, bottom=0.10, wspace=0.14)

        segment_boundaries = self._get_segment_boundaries(configs)
        last_cfg = configs[-1]

        show_anchors = any(cfg["show_anchors"] for cfg in configs)

        self._plot_slowness_track(
            axes[0],
            self.dtp,
            self.dtp_baseline,
            self.anchor_dtp,
            "P-wave slowness",
            "DTP (us/ft)",
            last_cfg["dtp_min"],
            last_cfg["dtp_max"],
            show_anchors,
            wave_key="DTP",
        )

        self._plot_slowness_track(
            axes[1],
            self.dts,
            self.dts_baseline,
            self.anchor_dts,
            "S-wave slowness",
            "DTS (us/ft)",
            last_cfg["dts_min"],
            last_cfg["dts_max"],
            show_anchors,
            wave_key="DTS",
        )

        self._plot_slowness_track(
            axes[2],
            self.dtst,
            self.dtst_baseline,
            self.anchor_dtst,
            "Stoneley-wave slowness",
            "DTST (us/ft)",
            last_cfg["dtst_min"],
            last_cfg["dtst_max"],
            show_anchors,
            wave_key="DTST",
        )

        self._plot_envelope_track(
            axes[3],
            self.s_dtp,
            "P-wave envelope",
            "S_DTP",
            last_cfg["s_dtp_threshold"],
            wave_key="DTP",
        )

        self._plot_envelope_track(
            axes[4],
            self.s_dts,
            "S-wave envelope",
            "S_DTS",
            last_cfg["s_dts_threshold"],
            wave_key="DTS",
        )

        self._plot_envelope_track(
            axes[5],
            self.s_dtst,
            "Stoneley-wave envelope",
            "S_DTST",
            last_cfg["s_dtst_threshold"],
            wave_key="DTST",
        )

        self._plot_lithology_track(axes[6])
        self._plot_profile_track(axes[7])

        for ax in axes:
            ax.set_ylim(np.min(self.depth), np.max(self.depth))
            ax.invert_yaxis()
            ax.tick_params(axis="y", labelsize=8)
            ax.tick_params(axis="x", labelsize=8, pad=3)
            ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=12))
            ax.grid(True, axis="y", alpha=0.18)

            for boundary in segment_boundaries:
                ax.axhline(
                    y=boundary,
                    color=(0 / 255, 0 / 255, 255 / 255),
                    linestyle="--",
                    linewidth=1.0,
                    alpha=0.85,
                )

        axes[0].set_ylabel("Depth (m)")

        self.embed_plot(fig)
        self.log_message("Integrated logging plot updated.", "SUCCESS")

    def _plot_slowness_track(
        self,
        ax,
        curve,
        baseline,
        anchor_mask,
        title,
        xlabel,
        x_min,
        x_max,
        show_anchors=False,
        wave_key="DTP",
    ):
        wave_color = WAVE_COLORS.get(wave_key, "black")
        baseline_color = BASELINE_COLORS.get(wave_key, "#1976d2")

        ax.plot(curve, self.depth, color=wave_color, linewidth=0.9, label="Original")
        ax.plot(
            baseline,
            self.depth,
            color=baseline_color,
            linewidth=0.85,
            linestyle="--",
            label="Floating baseline",
        )

        if show_anchors and anchor_mask is not None:
            valid_anchor = anchor_mask & np.isfinite(curve)
            ax.scatter(
                curve[valid_anchor],
                self.depth[valid_anchor],
                s=8,
                color="black",
                zorder=3,
                label="Anchors",
            )

        ax.set_title(title, fontsize=9)
        ax.set_xlabel(xlabel)

        if x_min != x_max:
            ax.set_xlim(x_max, x_min)

        ax.xaxis.set_ticks_position("bottom")
        ax.xaxis.set_label_position("bottom")
        ax.legend(fontsize=7, loc="best")

    def _plot_envelope_track(self, ax, envelope, title, xlabel, threshold, wave_key="DTP"):
        wave_color = WAVE_COLORS.get(wave_key, "#64b5f6")

        ax.plot(envelope, self.depth, color=wave_color, linewidth=0.9)

        if threshold > 0:
            ax.axvline(threshold, color="black", linestyle="--", linewidth=1.0)
            fill_where = envelope > threshold
            fill_start = threshold
        else:
            fill_where = np.isfinite(envelope)
            fill_start = 0.0

        ax.fill_betweenx(
            self.depth,
            fill_start,
            envelope,
            where=fill_where,
            color=wave_color,
            alpha=0.35,
        )

        valid = envelope[np.isfinite(envelope)]
        max_x = np.nanmax(valid) * 1.1 if valid.size > 0 else 1.0
        max_x = max(max_x, threshold * 1.2, 1.0)

        ax.set_title(title, fontsize=9)
        ax.set_xlabel(xlabel)
        ax.set_xlim(0, max_x)
        ax.xaxis.set_ticks_position("bottom")
        ax.xaxis.set_label_position("bottom")

    def _plot_lithology_track(self, ax):
        ax.set_title("FDA lithology", fontsize=9)
        ax.set_xlabel("Lithology")
        ax.set_xlim(0, 1)
        ax.set_xticks([])

        for i in range(len(self.depth) - 1):
            code = int(self.lithology[i])
            color = LITHOLOGY_COLORS.get(code, "white")

            ax.fill_betweenx(
                [self.depth[i], self.depth[i + 1]],
                0,
                1,
                color=color,
                alpha=0.95,
                linewidth=0,
            )

        labels = [
            "1 Sandstone",
            "2 Limestone",
            "3 Coal",
            "4 Mudstone",
        ]

        ax.text(
            0.5,
            0.02,
            "\n".join(labels),
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=7,
            bbox=dict(facecolor="white", alpha=0.72, edgecolor="none"),
        )

    def _plot_profile_track(self, ax):
        ax.set_title("Profile curves", fontsize=9)
        ax.set_xlabel("Content or value")

        if self.profiles is None or self.profiles.shape[1] == 0:
            ax.set_xlim(0, 1)
            ax.text(
                0.5,
                0.5,
                "No profile columns",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=8,
            )
            return

        cumulative = np.zeros_like(self.depth, dtype=float)

        for index in range(self.profiles.shape[1]):
            curve = self.profiles[:, index]
            next_curve = cumulative + curve
            color = PROFILE_DEFAULT_COLORS[index % len(PROFILE_DEFAULT_COLORS)]
            label = self.profile_names[index]

            ax.fill_betweenx(
                self.depth,
                cumulative,
                next_curve,
                where=np.isfinite(next_curve),
                color=color,
                alpha=0.90,
                label=label,
            )
            ax.plot(next_curve, self.depth, linewidth=0.35, color="black")
            cumulative = next_curve

        max_x = np.nanmax(cumulative) if np.any(np.isfinite(cumulative)) else 100.0
        max_x = max(1.0, max_x * 1.05)
        ax.set_xlim(0, max_x)
        ax.xaxis.set_ticks_position("bottom")
        ax.xaxis.set_label_position("bottom")
        ax.legend(fontsize=6, loc="best")

    def _get_segment_boundaries(self, configs):
        boundaries = set()

        for cfg in configs:
            top = min(cfg["top"], cfg["bottom"])
            bottom = max(cfg["top"], cfg["bottom"])
            boundaries.add(top)
            boundaries.add(bottom)

        boundaries = sorted(list(boundaries))

        if len(boundaries) >= 2:
            boundaries = boundaries[1:-1]

        return boundaries

    def embed_plot(self, fig):
        self.current_figure = fig

        layout = self.plot_frame.layout()
        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()

            if widget is not None:
                widget.setParent(None)

        container = QWidget()
        container_layout = QHBoxLayout(container)
        container_layout.setContentsMargins(0, 0, 0, 0)
        container_layout.setSpacing(4)

        plot_widget = QWidget()
        plot_layout = QVBoxLayout(plot_widget)
        plot_layout.setContentsMargins(0, 0, 0, 0)
        plot_layout.setSpacing(2)

        self.canvas = FigureCanvas(fig)
        self.toolbar = NavigationToolbar(self.canvas, self)

        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas)

        self.depth_scrollbar = QtWidgets.QScrollBar(QtCore.Qt.Vertical)
        self.depth_scrollbar.setFixedWidth(16)
        self.depth_scrollbar.valueChanged.connect(self.on_depth_scrollbar_changed)

        container_layout.addWidget(plot_widget, 1)
        container_layout.addWidget(self.depth_scrollbar)

        layout.addWidget(container)
        self.canvas.draw()

        self.full_depth_min = float(np.min(self.depth))
        self.full_depth_max = float(np.max(self.depth))
        self.current_view_min = self.full_depth_min
        self.current_view_max = self.full_depth_max

        self._update_depth_scrollbar()
        self.canvas.mpl_connect("scroll_event", self.on_canvas_scroll)

    def _set_depth_view(self, view_min, view_max):
        if self.depth is None or self.current_figure is None:
            return

        total_min = self.full_depth_min
        total_max = self.full_depth_max
        total_range = total_max - total_min

        if total_range <= 0:
            return

        min_view_range = max(total_range * 0.01, 0.5)

        if view_max < view_min:
            view_min, view_max = view_max, view_min

        view_range = view_max - view_min

        if view_range < min_view_range:
            center = (view_min + view_max) / 2
            view_min = center - min_view_range / 2
            view_max = center + min_view_range / 2
            view_range = min_view_range

        if view_range > total_range:
            view_min = total_min
            view_max = total_max
            view_range = total_range

        if view_min < total_min:
            view_min = total_min
            view_max = view_min + view_range

        if view_max > total_max:
            view_max = total_max
            view_min = view_max - view_range

        view_min = max(total_min, view_min)
        view_max = min(total_max, view_max)

        self.current_view_min = view_min
        self.current_view_max = view_max

        for ax in self.current_figure.axes:
            ax.set_ylim(view_max, view_min)

        self.canvas.draw_idle()
        self._update_depth_scrollbar()

    def _update_depth_scrollbar(self):
        if self.depth is None or self.current_figure is None or self.depth_scrollbar is None:
            return

        total_range = self.full_depth_max - self.full_depth_min
        view_range = self.current_view_max - self.current_view_min

        if total_range <= 0:
            self.depth_scrollbar.setRange(0, 0)
            return

        self.depth_scrollbar.blockSignals(True)

        if view_range >= total_range:
            self.depth_scrollbar.setRange(0, 0)
            self.depth_scrollbar.setPageStep(10000)
            self.depth_scrollbar.setValue(0)
        else:
            max_scroll = 10000
            self.depth_scrollbar.setRange(0, max_scroll)

            page_step = max(1, int(max_scroll * view_range / total_range))
            self.depth_scrollbar.setPageStep(page_step)

            denom = max(total_range - view_range, 1e-12)
            value = int(max_scroll * (self.current_view_min - self.full_depth_min) / denom)
            value = max(0, min(max_scroll, value))
            self.depth_scrollbar.setValue(value)

        self.depth_scrollbar.blockSignals(False)

    def on_depth_scrollbar_changed(self, value):
        if self.depth is None or self.depth_scrollbar is None:
            return

        total_min = self.full_depth_min
        total_range = self.full_depth_max - self.full_depth_min
        view_range = self.current_view_max - self.current_view_min

        if total_range <= 0 or view_range >= total_range:
            return

        max_scroll = self.depth_scrollbar.maximum()

        if max_scroll <= 0:
            return

        new_view_min = total_min + (total_range - view_range) * value / max_scroll
        new_view_max = new_view_min + view_range
        self._set_depth_view(new_view_min, new_view_max)

    def on_canvas_scroll(self, event):
        if self.current_figure is None or self.depth is None:
            return

        modifiers = QtWidgets.QApplication.keyboardModifiers()
        is_ctrl = bool(modifiers & QtCore.Qt.ControlModifier)

        view_min = self.current_view_min
        view_max = self.current_view_max
        view_range = view_max - view_min

        if event.ydata is not None:
            center = float(event.ydata)
        else:
            center = (view_min + view_max) / 2

        if is_ctrl:
            zoom_step = 0.12

            if event.button == "up":
                new_range = view_range * (1 - zoom_step)
            elif event.button == "down":
                new_range = view_range * (1 + zoom_step)
            else:
                return

            ratio = 0.5
            if view_range > 1e-12:
                ratio = (center - view_min) / view_range
                ratio = max(0.0, min(1.0, ratio))

            new_min = center - new_range * ratio
            new_max = new_min + new_range
            self._set_depth_view(new_min, new_max)
            return

        pan_step = view_range * 0.10

        if event.button == "up":
            new_min = view_min - pan_step
            new_max = view_max - pan_step
        elif event.button == "down":
            new_min = view_min + pan_step
            new_max = view_max + pan_step
        else:
            return

        self._set_depth_view(new_min, new_max)

    def out(self):
        if self.depth is None or self.s_dtp is None:
            self.log_message("No results are available. Please run the processing first.", "WARNING")
            return

        options = QFileDialog.Options()
        save_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Result File",
            "Slowness_Envelope_Result.txt",
            "Text Files (*.txt);;All Files (*)",
            options=options,
        )

        if not save_path:
            return

        try:
            columns = [
                self.depth,
                self.segment_id_result,
                self.dtp,
                self.dtp_baseline,
                self.s_dtp,
                self.dts,
                self.dts_baseline,
                self.s_dts,
                self.dtst,
                self.dtst_baseline,
                self.s_dtst,
                self.lithology,
            ]

            headers = [
                "Depth",
                "Segment_ID",
                "DTP",
                "DTP_BL",
                "S_DTP",
                "DTS",
                "DTS_BL",
                "S_DTS",
                "DTST",
                "DTST_BL",
                "S_DTST",
                "Lithology",
            ]

            if self.profiles is not None and self.profiles.shape[1] > 0:
                for index in range(self.profiles.shape[1]):
                    columns.append(self.profiles[:, index])
                    headers.append(self.profile_names[index])

            out_data = np.column_stack(columns)

            np.savetxt(
                save_path,
                out_data,
                fmt="%.8f",
                delimiter="\t",
                header="\t".join(headers),
                comments="",
            )

            self.log_message(f"Results exported successfully: {save_path}", "SUCCESS")

        except Exception as exc:
            self.log_message(f"Export failed: {exc}", "ERROR")
            QMessageBox.critical(self, "Export Error", str(exc))

    def save_figure(self):
        if self.current_figure is None:
            self.log_message("No figure is available. Please run the processing first.", "WARNING")
            return

        options = QFileDialog.Options()
        save_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Figure",
            "Slowness_Envelope_Figure.png",
            "PNG Files (*.png);;PDF Files (*.pdf);;TIFF Files (*.tif);;All Files (*)",
            options=options,
        )

        if not save_path:
            return

        try:
            self.current_figure.savefig(save_path, dpi=300, bbox_inches="tight")
            self.log_message(f"Figure exported successfully: {save_path}", "SUCCESS")
        except Exception as exc:
            self.log_message(f"Figure export failed: {exc}", "ERROR")
            QMessageBox.critical(self, "Figure Export Error", str(exc))


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())