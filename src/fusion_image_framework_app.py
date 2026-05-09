import sys
from pathlib import Path

import matplotlib
matplotlib.use("Qt5Agg")

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy.integrate import simpson

from PyQt5.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QSpinBox,
    QDoubleSpinBox,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas


DEFAULT_ATTRIBUTE_NAMES = ["S_DTP", "1-AMP", "S_DTS", "CAI", "S_DTST"]


def split_line(line: str) -> list[str]:
    line = line.strip()

    if "," in line:
        return [item.strip() for item in line.split(",") if item.strip()]

    if "\t" in line:
        return [item.strip() for item in line.split("\t") if item.strip()]

    return [item.strip() for item in line.split() if item.strip()]


def infer_delimiter(line: str):
    if "," in line:
        return ","
    if "\t" in line:
        return "\t"
    return None


def is_numeric_token(token: str) -> bool:
    try:
        float(token)
        return True
    except ValueError:
        return False


def read_table_with_optional_header(file_path: str | Path) -> tuple[np.ndarray, list[str]]:
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {file_path}")

    with file_path.open("r", encoding="utf-8-sig", errors="ignore") as file:
        lines = [line for line in file.readlines() if line.strip()]

    if not lines:
        raise ValueError("Input file is empty.")

    first_line = lines[0]
    first_tokens = split_line(first_line)
    delimiter = infer_delimiter(first_line)

    has_header = not all(is_numeric_token(token) for token in first_tokens)

    if has_header:
        column_names = first_tokens
        data = np.genfromtxt(file_path, delimiter=delimiter, skip_header=1)
    else:
        column_names = []
        data = np.genfromtxt(file_path, delimiter=delimiter, skip_header=0)

    if data.ndim == 1:
        data = data.reshape(1, -1)

    data = data[~np.isnan(data).all(axis=1)]

    if data.ndim != 2 or data.shape[0] == 0:
        raise ValueError("No valid numeric data were found in the input file.")

    if data.shape[1] < 2:
        raise ValueError("Input file must contain at least two columns: Depth and one attribute.")

    data = data[~np.isnan(data).any(axis=1)]

    if data.shape[0] < 3:
        raise ValueError("The input file contains too few valid numeric rows.")

    if has_header:
        if len(column_names) != data.shape[1]:
            raise ValueError(
                "The number of column names in the header does not match the number of numeric columns."
            )
    else:
        if data.shape[1] == 6:
            column_names = ["Depth"] + DEFAULT_ATTRIBUTE_NAMES
        else:
            column_names = ["Depth"] + [
                f"Attribute_{index + 1:02d}" for index in range(data.shape[1] - 1)
            ]

    return data, column_names


def load_attribute_data(file_path: str | Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    data, column_names = read_table_with_optional_header(file_path)

    depth = data[:, 0].astype(float)
    attributes = data[:, 1:].astype(float)
    attribute_names = column_names[1:]

    if np.any(~np.isfinite(depth)) or np.any(~np.isfinite(attributes)):
        raise ValueError("Input data contain NaN or infinite values.")

    if np.any(np.diff(depth) == 0):
        raise ValueError("Depth values contain duplicates.")

    if attributes.shape[1] < 1:
        raise ValueError("At least one attribute column is required.")

    return depth, attributes, attribute_names


def idw_interpolate_attributes(
    attributes: np.ndarray,
    n_interp: int = 50,
    power: float = 2.0,
) -> tuple[np.ndarray, np.ndarray]:
    if attributes.ndim != 2:
        raise ValueError("attributes must be a 2D array.")

    n_attributes = attributes.shape[1]

    if n_attributes < 1:
        raise ValueError("At least one attribute column is required.")

    if n_interp < n_attributes:
        raise ValueError(
            f"n_interp must be greater than or equal to the number of attributes ({n_attributes})."
        )

    if power <= 0:
        raise ValueError("power must be positive.")

    if n_attributes == 1:
        p_interp = np.linspace(1.0, 1.0, n_interp)
        fusion_image = np.repeat(attributes, n_interp, axis=1)
        return fusion_image, p_interp

    p_original = np.arange(1, n_attributes + 1, dtype=float)
    p_interp = np.linspace(1.0, float(n_attributes), n_interp)

    distance = np.abs(p_interp[None, :] - p_original[:, None])
    weights = np.zeros_like(distance, dtype=float)

    exact_mask = distance == 0
    non_exact_mask = ~exact_mask

    weights[non_exact_mask] = 1.0 / np.power(distance[non_exact_mask], power)

    for j in range(p_interp.size):
        if np.any(exact_mask[:, j]):
            weights[:, j] = 0.0
            exact_index = np.where(exact_mask[:, j])[0][0]
            weights[exact_index, j] = 1.0

    weights /= weights.sum(axis=0, keepdims=True)

    fusion_image = attributes @ weights
    return fusion_image, p_interp


def generate_fusion_image(
    attributes: np.ndarray,
    n_interp: int,
    idw_power: float,
    gaussian_sigma: float,
    clip_negative: bool,
) -> tuple[np.ndarray, np.ndarray]:
    image_interp, p_interp = idw_interpolate_attributes(
        attributes=attributes,
        n_interp=n_interp,
        power=idw_power,
    )

    image_smooth = gaussian_filter(image_interp, sigma=gaussian_sigma)

    if clip_negative:
        image_smooth = np.maximum(image_smooth, 0.0)

    return image_smooth, p_interp


def compute_composite_indicator(
    depth: np.ndarray,
    fusion_image: np.ndarray,
    p_interp: np.ndarray,
    window_size: int,
    step: int,
    clip_negative: bool,
) -> tuple[np.ndarray, np.ndarray]:
    if window_size < 2:
        raise ValueError("window_size must be greater than or equal to 2.")

    if step < 1:
        raise ValueError("step must be greater than or equal to 1.")

    if len(depth) != fusion_image.shape[0]:
        raise ValueError("depth length must match the number of rows in fusion_image.")

    if window_size > len(depth):
        raise ValueError("window_size cannot exceed the number of depth samples.")

    indicator_values = []
    indicator_depths = []

    for start in range(0, len(depth) - window_size + 1, step):
        end = start + window_size

        window_image = fusion_image[start:end, :]
        window_depth = depth[start:end]

        area_along_attributes = simpson(window_image, x=p_interp, axis=1)
        volume = simpson(area_along_attributes, x=window_depth)

        center_index = start + window_size // 2

        indicator_values.append(volume)
        indicator_depths.append(depth[center_index])

    indicator_values = np.asarray(indicator_values, dtype=float)
    indicator_depths = np.asarray(indicator_depths, dtype=float)

    if clip_negative:
        indicator_values = np.maximum(indicator_values, 0.0)

    return indicator_depths, indicator_values


def build_full_length_indicator(
    depth: np.ndarray,
    indicator_depths: np.ndarray,
    indicator_values: np.ndarray,
) -> np.ndarray:
    full_indicator = np.full(depth.shape, np.nan, dtype=float)

    for indicator_depth, indicator_value in zip(indicator_depths, indicator_values):
        nearest_index = np.argmin(np.abs(depth - indicator_depth))
        full_indicator[nearest_index] = indicator_value

    return full_indicator


def export_results(
    output_path: str | Path,
    depth: np.ndarray,
    full_indicator: np.ndarray,
    fusion_image: np.ndarray,
) -> None:
    output_path = Path(output_path)

    output_data = np.column_stack((depth, full_indicator, fusion_image))
    output_data[:, 1:] = np.maximum(output_data[:, 1:], 0.0)

    headers = ["Depth", "LF"] + [
        f"Img_{index + 1:02d}" for index in range(fusion_image.shape[1])
    ]

    np.savetxt(
        output_path,
        output_data,
        fmt="%.6f",
        delimiter="\t",
        header="\t".join(headers),
        comments="",
    )


class FusionImageCanvas(FigureCanvas):
    def __init__(self, parent=None):
        self.fig, self.axes = plt.subplots(1, 2, figsize=(12, 8), sharey=True)
        super().__init__(self.fig)
        self.setParent(parent)
        self.colorbar = None

    def clear(self):
        self.fig.clear()
        self.axes = self.fig.subplots(1, 2, sharey=True)
        self.colorbar = None
        self.draw()

    def plot_results(
        self,
        depth: np.ndarray,
        fusion_image: np.ndarray,
        p_interp: np.ndarray,
        indicator_depths: np.ndarray,
        indicator_values: np.ndarray,
        cmap: str,
    ):
        self.fig.clear()
        axes = self.fig.subplots(1, 2, sharey=True)

        image = axes[0].imshow(
            fusion_image,
            aspect="auto",
            extent=[p_interp[0], p_interp[-1], depth[-1], depth[0]],
            origin="upper",
            cmap=cmap,
        )

        axes[0].set_xlabel("Attribute axis")
        axes[0].set_ylabel("Depth")
        axes[0].set_title("Multi-attribute fusion image")

        axes[1].plot(
            indicator_values,
            indicator_depths,
            color="black",
            linewidth=1.4,
        )

        axes[1].set_xlabel("Composite fracture indicator")
        axes[1].set_title("Composite fracture indicator curve")
        axes[1].set_xlim(left=0)
        axes[1].grid(True, linewidth=0.5, alpha=0.6)
        axes[1].invert_yaxis()

        colorbar_axis = self.fig.add_axes([0.08, 0.93, 0.35, 0.025])
        self.fig.colorbar(image, cax=colorbar_axis, orientation="horizontal")

        self.fig.tight_layout(rect=[0, 0, 1, 0.9])
        self.draw()


class FusionImageApp(QWidget):
    def __init__(self):
        super().__init__()

        self.input_path: Path | None = None
        self.depth: np.ndarray | None = None
        self.attributes: np.ndarray | None = None
        self.attribute_names: list[str] = []
        self.fusion_image: np.ndarray | None = None
        self.p_interp: np.ndarray | None = None
        self.indicator_depths: np.ndarray | None = None
        self.indicator_values: np.ndarray | None = None
        self.full_indicator: np.ndarray | None = None

        self.init_ui()

    def init_ui(self):
        self.setWindowTitle("Array Acoustic Multi-attribute Fusion Image")
        self.resize(1350, 850)

        main_layout = QHBoxLayout(self)

        left_panel = QVBoxLayout()
        left_panel.setSpacing(10)

        self.file_line = QLineEdit()
        self.file_line.setReadOnly(True)

        browse_button = QPushButton("Select Input File")
        browse_button.clicked.connect(self.select_file)

        file_group = QGroupBox("Input")
        file_layout = QVBoxLayout()
        file_layout.addWidget(QLabel("Input column order:"))
        file_layout.addWidget(QLabel("Depth, Attribute_1, Attribute_2, ..."))
        file_layout.addWidget(QLabel("If a header row is provided, attribute names are read automatically."))
        file_layout.addWidget(self.file_line)
        file_layout.addWidget(browse_button)
        file_group.setLayout(file_layout)

        parameter_group = QGroupBox("Parameters")
        parameter_layout = QFormLayout()

        self.n_interp_spin = QSpinBox()
        self.n_interp_spin.setRange(1, 2000)
        self.n_interp_spin.setValue(50)

        self.idw_power_spin = QDoubleSpinBox()
        self.idw_power_spin.setRange(0.1, 10.0)
        self.idw_power_spin.setSingleStep(0.1)
        self.idw_power_spin.setValue(2.0)

        self.gaussian_sigma_spin = QDoubleSpinBox()
        self.gaussian_sigma_spin.setRange(0.0, 20.0)
        self.gaussian_sigma_spin.setSingleStep(0.1)
        self.gaussian_sigma_spin.setValue(1.0)

        self.window_size_spin = QSpinBox()
        self.window_size_spin.setRange(2, 10000)
        self.window_size_spin.setValue(10)

        self.step_spin = QSpinBox()
        self.step_spin.setRange(1, 1000)
        self.step_spin.setValue(1)

        self.clip_check = QCheckBox("Clip negative values to zero")
        self.clip_check.setChecked(True)

        self.cmap_combo = QComboBox()
        self.cmap_combo.addItems(["jet", "viridis", "turbo", "plasma", "inferno", "magma"])
        self.cmap_combo.setCurrentText("jet")

        parameter_layout.addRow("Interpolated columns", self.n_interp_spin)
        parameter_layout.addRow("IDW power", self.idw_power_spin)
        parameter_layout.addRow("Gaussian sigma", self.gaussian_sigma_spin)
        parameter_layout.addRow("Window size", self.window_size_spin)
        parameter_layout.addRow("Step", self.step_spin)
        parameter_layout.addRow("Colormap", self.cmap_combo)
        parameter_layout.addRow(self.clip_check)

        parameter_group.setLayout(parameter_layout)

        run_button = QPushButton("Run")
        run_button.clicked.connect(self.run_processing)
        run_button.setMinimumHeight(36)

        export_button = QPushButton("Export Results")
        export_button.clicked.connect(self.export_result_file)
        export_button.setMinimumHeight(32)

        save_fig_button = QPushButton("Save Figure")
        save_fig_button.clicked.connect(self.save_figure)
        save_fig_button.setMinimumHeight(32)

        self.log_box = QTextEdit()
        self.log_box.setReadOnly(True)
        self.log_box.setMinimumHeight(220)

        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)

        left_panel.addWidget(file_group)
        left_panel.addWidget(parameter_group)
        left_panel.addWidget(run_button)
        left_panel.addWidget(export_button)
        left_panel.addWidget(save_fig_button)
        left_panel.addWidget(line)
        left_panel.addWidget(QLabel("Log"))
        left_panel.addWidget(self.log_box)
        left_panel.addStretch()

        self.canvas = FusionImageCanvas(self)

        main_layout.addLayout(left_panel, 0)
        main_layout.addWidget(self.canvas, 1)

    def log(self, message: str):
        self.log_box.append(message)

    def select_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select input attribute file",
            "",
            "Data files (*.txt *.csv *.dat);;Text files (*.txt);;CSV files (*.csv);;All files (*.*)",
        )

        if not file_path:
            return

        self.input_path = Path(file_path)
        self.file_line.setText(str(self.input_path))

        try:
            self.depth, self.attributes, self.attribute_names = load_attribute_data(self.input_path)

            if self.n_interp_spin.value() < self.attributes.shape[1]:
                self.n_interp_spin.setValue(self.attributes.shape[1])

            self.log(f"Loaded file: {self.input_path}")
            self.log(f"Number of samples: {len(self.depth)}")
            self.log(f"Number of attributes: {self.attributes.shape[1]}")
            self.log(f"Attribute order: {', '.join(self.attribute_names)}")
        except Exception as exc:
            self.depth = None
            self.attributes = None
            self.attribute_names = []
            QMessageBox.critical(self, "Load Error", str(exc))
            self.log(f"Load error: {exc}")

    def run_processing(self):
        if self.input_path is None:
            QMessageBox.warning(self, "No Input File", "Please select an input file first.")
            return

        if self.depth is None or self.attributes is None:
            try:
                self.depth, self.attributes, self.attribute_names = load_attribute_data(self.input_path)
            except Exception as exc:
                QMessageBox.critical(self, "Load Error", str(exc))
                return

        try:
            n_interp = self.n_interp_spin.value()
            idw_power = self.idw_power_spin.value()
            gaussian_sigma = self.gaussian_sigma_spin.value()
            window_size = self.window_size_spin.value()
            step = self.step_spin.value()
            clip_negative = self.clip_check.isChecked()

            if n_interp < self.attributes.shape[1]:
                raise ValueError(
                    f"Interpolated columns ({n_interp}) must be greater than or equal to "
                    f"the number of attributes ({self.attributes.shape[1]})."
                )

            if window_size > len(self.depth):
                raise ValueError(
                    f"Window size ({window_size}) cannot exceed the number of depth samples ({len(self.depth)})."
                )

            self.fusion_image, self.p_interp = generate_fusion_image(
                attributes=self.attributes,
                n_interp=n_interp,
                idw_power=idw_power,
                gaussian_sigma=gaussian_sigma,
                clip_negative=clip_negative,
            )

            self.indicator_depths, self.indicator_values = compute_composite_indicator(
                depth=self.depth,
                fusion_image=self.fusion_image,
                p_interp=self.p_interp,
                window_size=window_size,
                step=step,
                clip_negative=clip_negative,
            )

            self.full_indicator = build_full_length_indicator(
                depth=self.depth,
                indicator_depths=self.indicator_depths,
                indicator_values=self.indicator_values,
            )

            self.canvas.plot_results(
                depth=self.depth,
                fusion_image=self.fusion_image,
                p_interp=self.p_interp,
                indicator_depths=self.indicator_depths,
                indicator_values=self.indicator_values,
                cmap=self.cmap_combo.currentText(),
            )

            self.log("Processing completed.")
            self.log(f"Number of attributes: {self.attributes.shape[1]}")
            self.log(f"Interpolated columns: {n_interp}")
            self.log(f"IDW power: {idw_power}")
            self.log(f"Gaussian sigma: {gaussian_sigma}")
            self.log(f"Window size: {window_size}")
            self.log(f"Step: {step}")

        except Exception as exc:
            QMessageBox.critical(self, "Processing Error", str(exc))
            self.log(f"Processing error: {exc}")

    def export_result_file(self):
        if self.depth is None or self.full_indicator is None or self.fusion_image is None:
            QMessageBox.warning(self, "No Results", "Please run the processing first.")
            return

        default_path = self.input_path.with_name(
            f"{self.input_path.stem}_fusion_image_output.txt"
        )

        output_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save output data",
            str(default_path),
            "Text files (*.txt);;CSV files (*.csv);;All files (*.*)",
        )

        if not output_path:
            return

        try:
            export_results(
                output_path=output_path,
                depth=self.depth,
                full_indicator=self.full_indicator,
                fusion_image=self.fusion_image,
            )
            self.log(f"Results exported to: {output_path}")
            QMessageBox.information(self, "Export Complete", "Results exported successfully.")
        except Exception as exc:
            QMessageBox.critical(self, "Export Error", str(exc))
            self.log(f"Export error: {exc}")

    def save_figure(self):
        if self.fusion_image is None:
            QMessageBox.warning(self, "No Figure", "Please run the processing first.")
            return

        default_path = self.input_path.with_name(
            f"{self.input_path.stem}_fusion_image_figure.png"
        )

        figure_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save figure",
            str(default_path),
            "PNG files (*.png);;PDF files (*.pdf);;TIFF files (*.tif);;All files (*.*)",
        )

        if not figure_path:
            return

        try:
            self.canvas.fig.savefig(figure_path, dpi=300, bbox_inches="tight")
            self.log(f"Figure saved to: {figure_path}")
            QMessageBox.information(self, "Save Complete", "Figure saved successfully.")
        except Exception as exc:
            QMessageBox.critical(self, "Save Error", str(exc))
            self.log(f"Figure save error: {exc}")


def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")

    window = FusionImageApp()
    window.show()

    sys.exit(app.exec_())


if __name__ == "__main__":
    main()