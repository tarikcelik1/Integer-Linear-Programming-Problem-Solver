import sys
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, QWidget, QLabel, QPushButton, QLineEdit, 
    QComboBox, QScrollArea, QMessageBox, QTableWidget, QTableWidgetItem, QHeaderView, QGroupBox
)
from PyQt6.QtCore import Qt

from branch_and_bound_max import simplex, branch_and_bound
from gomory import GomoryCut

class ILPSolverApp(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Integer Linear Programming Solver")
        self.setGeometry(100, 100, 1000, 600)

        self.main_widget = QWidget()
        self.setCentralWidget(self.main_widget)

        self.layout = QVBoxLayout()
        self.main_widget.setLayout(self.layout)

        self.header = QLabel("Linear Programming Solver", self)
        self.header.setStyleSheet("font-size: 20px; font-weight: bold; text-align: center;")
        self.header.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.layout.addWidget(self.header)

        self.input_section = QWidget()
        self.input_layout = QVBoxLayout()
        self.input_section.setLayout(self.input_layout)
        self.layout.addWidget(self.input_section)

        self.results_section = QScrollArea()
        self.results_section.setWidgetResizable(True)
        self.results_widget = QWidget()
        self.results_layout = QVBoxLayout()
        self.results_widget.setLayout(self.results_layout)
        self.results_section.setWidget(self.results_widget)
        self.layout.addWidget(self.results_section)

        self.setup_inputs()
        self.show_help_button()

    def setup_inputs(self):
        # Number of Variables and Constraints
        input_group = QGroupBox("Setup Problem")
        input_group_layout = QHBoxLayout()
        input_group.setLayout(input_group_layout)

        self.num_vars_label = QLabel("Number of Variables:")
        self.num_vars_input = QLineEdit()
        self.num_vars_input.setPlaceholderText("Enter number of variables")
        self.num_vars_input.setFixedWidth(80)

        self.num_constraints_label = QLabel("Number of Constraints:")
        self.num_constraints_input = QLineEdit()
        self.num_constraints_input.setPlaceholderText("Enter number of constraints")
        self.num_constraints_input.setFixedWidth(80)

        self.next_button = QPushButton("Next")
        self.next_button.setStyleSheet("font-weight: bold; padding: 8px; font-size: 12px;")
        self.next_button.clicked.connect(self.create_input_fields)

        input_group_layout.addWidget(self.num_vars_label)
        input_group_layout.addWidget(self.num_vars_input)
        input_group_layout.addWidget(self.num_constraints_label)
        input_group_layout.addWidget(self.num_constraints_input)
        input_group_layout.addWidget(self.next_button)

        self.input_layout.addWidget(input_group)

    def show_help_button(self):
        help_button = QPushButton("Help")
        help_button.setStyleSheet("font-weight: bold; padding: 8px; font-size: 12px;")
        help_button.clicked.connect(self.show_help)
        self.input_layout.addWidget(help_button)

    def create_input_fields(self):
        try:
            num_vars = int(self.num_vars_input.text())
            num_constraints = int(self.num_constraints_input.text())

            if num_vars <= 0 or num_constraints <= 0:
                raise ValueError("Number of variables and constraints must be positive integers.")

            # Clear Results Section
            for i in reversed(range(self.results_layout.count())):
                self.results_layout.itemAt(i).widget().deleteLater()

            # Objective Function Input
            obj_group = QGroupBox("Objective Function")
            obj_layout = QHBoxLayout()
            obj_group.setLayout(obj_layout)
            self.results_layout.addWidget(obj_group)

            self.obj_inputs = []
            for i in range(num_vars):
                label = QLabel(f"x{i+1}:")
                label.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
                input_field = QLineEdit()
                input_field.setFixedWidth(50)

                obj_layout.addWidget(label)
                obj_layout.addWidget(input_field)
                self.obj_inputs.append(input_field)

            # Constraints Input
            self.constraints_inputs = []

            for i in range(num_constraints):
                constraint_group = QGroupBox(f"Constraint {i + 1}")
                constraint_layout = QHBoxLayout()
                constraint_group.setLayout(constraint_layout)
                self.results_layout.addWidget(constraint_group)

                row_inputs = []
                for j in range(num_vars):
                    label = QLabel(f"x{j+1}:")
                    label.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
                    input_field = QLineEdit()
                    input_field.setFixedWidth(50)

                    constraint_layout.addWidget(label)
                    constraint_layout.addWidget(input_field)
                    row_inputs.append(input_field)

                sense_selector = QComboBox()
                sense_selector.addItems(["<=", ">=", "="])
                rhs_input = QLineEdit()
                rhs_input.setFixedWidth(80)

                constraint_layout.addWidget(sense_selector)
                constraint_layout.addWidget(QLabel("RHS:"))
                constraint_layout.addWidget(rhs_input)

                self.constraints_inputs.append((row_inputs, sense_selector, rhs_input))

            # Solve Button
            solve_button = QPushButton("Solve")
            solve_button.setStyleSheet("font-weight: bold; padding: 8px; font-size: 12px;")
            solve_button.clicked.connect(self.process_inputs)
            self.results_layout.addWidget(solve_button)

        except ValueError as e:
            QMessageBox.warning(self, "Input Error", str(e))

    def process_inputs(self):
        try:
            c = [float(input_field.text()) for input_field in self.obj_inputs]

            A = []
            b = []
            constraint_types = []

            for row_inputs, sense_selector, rhs_input in self.constraints_inputs:
                row = [float(input_field.text()) for input_field in row_inputs]
                sense = sense_selector.currentText()
                rhs = float(rhs_input.text())

                A.append(row)
                b.append(rhs)
                constraint_types.append(sense)

            # Example solving process (replace with actual solver logic)
            solutions, obj_values = simplex(c, A, b)
            simplex_solution = f"Simplex Solution: {solutions} Objective Value: {obj_values}"

            solutionbb, obj_valuebb = branch_and_bound(c, A, b, constraint_types)
            branch_and_bound_solution = f"Branch and Bound Solution: {solutionbb} Objective Value: {obj_valuebb}"

            gc = GomoryCut(A, b, c, constraint_types)
            integer_solutiongc, integer_obj_valuegc = gc.solve()
            gomory_cut_solution = f"Gomory Cut Solution: {integer_solutiongc} Objective Value: {integer_obj_valuegc}"

            # Display Results
            self.show_results(simplex_solution, branch_and_bound_solution, gomory_cut_solution)

        except ValueError as e:
            QMessageBox.warning(self, "Input Error", str(e))

    def show_results(self, lp_solution, ilp_solution, ilp_solutiongc):
        for i in reversed(range(self.results_layout.count())):
            self.results_layout.itemAt(i).widget().deleteLater()

        results_table = QTableWidget(3, 2)
        results_table.setHorizontalHeaderLabels(["Method", "Solution"])
        results_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        results_table.setItem(0, 0, QTableWidgetItem("Simplex"))
        results_table.setItem(0, 1, QTableWidgetItem(lp_solution))
        results_table.setItem(1, 0, QTableWidgetItem("Branch and Bound"))
        results_table.setItem(1, 1, QTableWidgetItem(ilp_solution))
        results_table.setItem(2, 0, QTableWidgetItem("Gomory Cut"))
        results_table.setItem(2, 1, QTableWidgetItem(ilp_solutiongc))

        self.results_layout.addWidget(results_table)

    def show_help(self):
        help_text = (
            "Input Format:\n\n"
            "Objective Function:\n"
            "Enter coefficients of the variables (e.g., if the objective function is 3x1 + 5x2, enter 3 5).\n\n"
            "Constraints:\n"
            "For each constraint, specify the coefficients of the variables, the relation (<=, >=, =), "
            "and the RHS value.\n\n"
            "Example:\n"
            "Objective Function: Maximize 3x1 + 5x2\n"
            "Constraints:\n"
            " 1. 2x1 + 3x2 <= 10\n"
            " 2. x1 + x2 >= 4\n"
            " 3. x1 = 2\n"
        )
        QMessageBox.information(self, "Input Help", help_text)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = ILPSolverApp()
    window.show()
    sys.exit(app.exec())
