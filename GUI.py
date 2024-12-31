import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from functools import partial
from branch_and_bound_max import simplex, branch_and_bound
from gomory import GomoryCut

def show_help():
    help_text = (
        "Input Format:\n\n"
        "Objective Function:\n"
        "Enter coefficients of the variables (e.g., if the objective function is 3x1 + 5x2, enter 3 and 5).\n\n"
        "Constraints:\n"
        "For each constraint, specify the coefficients of the variables, the relation (<=, >=, =), "
        "and the RHS value (e.g., 2x1 + 3x2 <= 10).\n\n"
        "Example:\n"
        "Objective Function: Maximize 3x1 + 5x2\n"
        "Constraints:\n"
        " 1. 2x1 + 3x2 <= 10\n"
        " 2. x1 + x2 >= 4\n"
        " 3. x1 = 2\n"
    )
    messagebox.showinfo("Input Help", help_text)

def show_error(msg):
    # Instead of messageboxes, we show errors in the results frame
    for widget in results_frame.winfo_children():
        widget.destroy()
    ttk.Label(results_frame, text="Error:", font=('Arial', 12, 'bold'), foreground="red").grid(row=0, column=0, sticky="w")
    ttk.Label(results_frame, text=msg, font=('Arial', 11)).grid(row=1, column=0, sticky="w")

def show_results(lp_solution, ilp_solution, ilp_solutiongc):
    # Clear previous results
    for widget in results_frame.winfo_children():
        widget.destroy()
    ttk.Label(results_frame, text="Results", font=('Arial', 14, 'bold')).grid(row=0, column=0, sticky="w", pady=(10, 5))
    ttk.Label(results_frame, text=lp_solution, font=('Arial', 12)).grid(row=1, column=0, sticky="w", padx=10, pady=5)
    ttk.Label(results_frame, text=ilp_solution, font=('Arial', 12)).grid(row=2, column=0, sticky="w", padx=10, pady=5)
    ttk.Label(results_frame, text=ilp_solutiongc, font=('Arial', 12)).grid(row=3, column=0, sticky="w", padx=10, pady=5)

def solve_problem():
    try:
        num_vars = int(entry_vars.get())
        num_constraints = int(entry_constraints.get())
        
        if num_vars <= 0 or num_constraints <= 0:
            raise ValueError("Number of variables and constraints must be positive integers.")
        
        # Clear the input frame for new fields
        for widget in input_frame.winfo_children():
            widget.destroy()

        ttk.Label(input_frame, text="Enter Coefficients for Objective Function:",
                  font=('Arial', 12, 'bold')).grid(row=1, column=0, columnspan=2, pady=(10, 10), sticky="w")

        obj_entries = []
        obj_frame = ttk.Frame(input_frame)
        obj_frame.grid(row=2, column=0, columnspan=2, pady=5, sticky="w")

        for i in range(num_vars):
            ttk.Label(obj_frame, text=f"x{i+1}:", font=('Arial', 10)).grid(row=0, column=i*2, padx=5, pady=5, sticky="e")
            oe = ttk.Entry(obj_frame, width=5)
            oe.grid(row=0, column=i*2+1, padx=5, pady=5)
            obj_entries.append(oe)

        # Constraints Section
        ttk.Label(input_frame, text="Constraints:", font=('Arial', 12, 'bold')).grid(row=3, column=0, columnspan=2, pady=(20, 10), sticky="w")

        constraint_entries = []
        senses = []
        rhs_entries = []

        for i in range(num_constraints):
            frame_label = ttk.LabelFrame(input_frame, text=f"Constraint {i+1}", padding=10)
            frame_label.grid(row=4+i, column=0, columnspan=2, pady=10, sticky="w")

            row_entries = []
            col_index = 0
            for j in range(num_vars):
                ttk.Label(frame_label, text=f"x{j+1}:", font=('Arial', 10)).grid(row=0, column=col_index, padx=5, pady=5, sticky="e")
                coef_entry = ttk.Entry(frame_label, width=5)
                coef_entry.grid(row=0, column=col_index+1, padx=5, pady=5)
                row_entries.append(coef_entry)
                col_index += 2

            sense_var = tk.StringVar(value="<=")
            sense_combo = ttk.Combobox(frame_label, textvariable=sense_var, values=["<=", ">=", "="], width=5, state="readonly")
            sense_combo.grid(row=0, column=col_index, padx=5, pady=5)
            senses.append(sense_var)

            ttk.Label(frame_label, text="RHS:", font=('Arial', 10)).grid(row=0, column=col_index+1, padx=5, pady=5, sticky="e")
            rhs_entry = ttk.Entry(frame_label, width=5)
            rhs_entry.grid(row=0, column=col_index+2, padx=5, pady=5)
            rhs_entries.append(rhs_entry)

            constraint_entries.append(row_entries)

        solve_btn = ttk.Button(input_frame, text="Solve", 
                               command=partial(process_inputs, num_vars, num_constraints, obj_entries, constraint_entries, senses, rhs_entries, problem_type_var))
        solve_btn.grid(row=4+num_constraints, column=0, columnspan=2, pady=20)

    except ValueError as e:
        show_error(str(e))


def process_inputs(num_vars, num_constraints, obj_entries, constraint_entries, senses, rhs_entries, problem_type_var):
    try:
        problem_type = problem_type_var.get()
        c = []
        for i, obj_entry in enumerate(obj_entries):
            val = obj_entry.get().strip()
            if not val:
                raise ValueError(f"Objective function coefficient for x{i+1} cannot be empty.")
            c.append(float(val))

        A = []
        b = []
        constraint_type = [None] * num_constraints
        for i in range(num_constraints):
            row_coefs = []
            for entry in constraint_entries[i]:
                val = entry.get().strip()
                if not val:
                    raise ValueError(f"Constraint coefficient for constraint {i+1} cannot be empty.")
                row_coefs.append(float(val))

            rhs_val = rhs_entries[i].get().strip()
            if not rhs_val:
                raise ValueError(f"RHS for constraint {i+1} cannot be empty.")
            rhs_num = float(rhs_val)
            sense = senses[i].get()
            constraint_type[i] = sense
            if sense == ">=":
                row_coefs = [coef for coef in row_coefs]
                rhs_num = rhs_num
            elif sense == "=":
                pass
            A.append(row_coefs)
            b.append(rhs_num)
        print("A:",A)
        print("b:",b)

        #Multiply c by -1 for simplex
        print(c)
        solution, obj_value = simplex(c, A, b)
        

        integer_solutionbb, integer_obj_valuebb = branch_and_bound(c, A, b,constraint_type)

        gc = GomoryCut(A, b, c,constraint_type)
        integer_solutiongc, integer_obj_valuegc = gc.solve()
        #take the solution with 2 digits after the decimal point
        #convert solution to integer values
        integer_solutiongc = [round(x, 2) for x in integer_solutiongc]
        integer_solutiongc = [float(x) for x in integer_solutiongc]
        integer_obj_valuegc = round(integer_obj_valuegc, 2)

        lpsol = ""
        ilpsol = ""
        ilpsolgc = ""
        for i in range(len(solution)):
            if i !=0:
                lpsol += ", "
            lpsol += "x"+str(i+1)+"="+str(solution[i])
        for i in range(len(integer_solutionbb)):
            if i !=0:
                ilpsol += ", "
            ilpsol += "x"+str(i+1)+"="+str(integer_solutionbb[i])
        for i in range(len(integer_solutiongc)):
            if i !=0:
                ilpsolgc += ", "
            ilpsolgc += "x"+str(i+1)+"="+str(integer_solutiongc[i])

        lp_solution = f"Simplex Solution: {lpsol}, Optimal Value: {obj_value}"
        ilp_solutionbb = f"Branch and Bound Solution: {ilpsol}, Optimal Value: {integer_obj_valuebb}"
        ilp_solutiongc = f"Gomory Cut Solution: {ilpsolgc}, Optimal Value: {integer_obj_valuegc}"

        show_results(lp_solution, ilp_solutionbb, ilp_solutiongc)

    except ValueError as e:
        show_error(f"Invalid input: {e}")
    except Exception as e:
        show_error(f"Error: {e}")


# GUI Setup
root = tk.Tk()
root.title("Linear Programming Solver")

problem_type_var = tk.StringVar()

style = ttk.Style(root)
style.theme_use('default')
style.configure('.', font=('Arial', 10))
style.configure('TLabel', foreground='#333')
style.configure('TEntry', padding=5)
style.configure('TButton', padding=5)

main_frame = ttk.Frame(root, padding=20)
main_frame.pack(fill="both", expand=True)

top_frame = ttk.Frame(main_frame)
top_frame.pack(side="top", fill="x")

ttk.Button(top_frame, text="Help", command=show_help).grid(row=0, column=5, padx=10, pady=5)

ttk.Label(top_frame, text="Number of Variables:", font=('Arial', 10, 'bold')).grid(row=0, column=0, padx=5, pady=5, sticky="e")
entry_vars = ttk.Entry(top_frame, width=5)
entry_vars.grid(row=0, column=1, padx=5, pady=5, sticky="w")

ttk.Label(top_frame, text="Number of Constraints:", font=('Arial', 10, 'bold')).grid(row=0, column=2, padx=5, pady=5, sticky="e")
entry_constraints = ttk.Entry(top_frame, width=5)
entry_constraints.grid(row=0, column=3, padx=5, pady=5, sticky="w")

ttk.Button(top_frame, text="Next", command=solve_problem).grid(row=0, column=4, padx=10, pady=5)

input_frame = ttk.Frame(main_frame)
input_frame.pack(side="top", fill="both", expand=True, pady=10)

results_frame = ttk.Frame(main_frame)
results_frame.pack(side="bottom", fill="x", pady=10)

root.mainloop()
