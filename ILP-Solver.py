import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from functools import partial
from branch_and_bound import simplex, branch_and_bound
from gomory import GomoryCut


font_type = "Helvetica" # Font type for the GUI


def to_subscript(num):
    subscript_map = {
        '0': '₀', '1': '₁', '2': '₂', '3': '₃', '4': '₄',
        '5': '₅', '6': '₆', '7': '₇', '8': '₈', '9': '₉'
    }
    return ''.join(subscript_map[digit] for digit in str(num))

def show_help():
    help_text = (
        "Input Format:\n\n"
        "Objective Function:\n"
        "Enter coefficients of the variables (e.g., if the objective function is 3x1 + 5x2, enter 3 and 5).\n\n"
        "Constraints:\n"
        "and the RHS value (e.g., 2x1 + 3x2 <= 10).\n\n"
        "Example:\n"
        "Objective Function: Maximize 3x1 + 5x2\n"
        "Constraints:\n"
        " 1. 2x1 + 3x2 <= 10\n"
        " 2. x1 + x2 <= 4\n"
        " 3. x1 <= 2\n"
    )
    messagebox.showinfo("Input Help", help_text)

def show_error(msg):
    # Instead of messageboxes, we show errors in the results frame
    for widget in results_frame.winfo_children():
        widget.destroy()
    ttk.Label(results_frame, text="Error:", font=(font_type, 12, 'bold'), foreground="red").grid(row=0, column=0, sticky="w")
    ttk.Label(results_frame, text=msg, font=(font_type, 11)).grid(row=1, column=0, sticky="w")

def show_results(lp_solution, ilp_solution, ilp_solutiongc):
    # Clear previous results
    for widget in results_frame.winfo_children():
        widget.destroy()
    tk.Label(results_frame, text="Results",bg=root_color, font=(font_type, 16, 'bold')).grid(row=0, column=0, sticky="w", pady=(10, 5))
    #ttk.Label(results_frame, text=lp_solution, font=(font_type, 12)).grid(row=1, column=0, sticky="w", padx=10, pady=5)
    tk.Label(results_frame, text=ilp_solution, bg=root_color, font=(font_type, 14)).grid(row=1, column=0, sticky="w", padx=10, pady=5)
    tk.Label(results_frame, text=ilp_solutiongc, bg=root_color, font=(font_type, 14)).grid(row=2, column=0, sticky="w", padx=10, pady=5)

def solve_problem():
    try:
        num_vars = int(entry_vars.get())
        num_constraints = int(entry_constraints.get())
        
        if num_vars <= 0 or num_constraints <= 0:
            raise ValueError("Number of variables and constraints must be positive integers.")
        
        # Clear the input frame for new fields
        for widget in input_frame.winfo_children():
            widget.destroy()
        tk.Label(input_frame, text="Enter Coefficients for Objective Function:",bg=root_color,
                  font=(font_type, 14, 'bold')).grid(row=1, column=0, columnspan=2, pady=(10, 10), sticky="w")

        obj_entries = []
        obj_frame = ttk.Frame(input_frame)
        obj_frame.grid(row=2, column=0, columnspan=2, pady=5, sticky="w")

        for i in range(num_vars):
            ttk.Label(obj_frame, text=f"x" + to_subscript(i + 1), font=(font_type, 13)).grid(row=0, column=i*2, padx=5, pady=5, sticky="e")
            oe = ttk.Entry(obj_frame, width=10)
            oe.grid(row=0, column=i*2+1, padx=5, pady=5)
            obj_entries.append(oe)

        # Constraints Section
        tk.Label(input_frame, text="Constraints:", bg=root_color, font=(font_type, 14, 'bold')).grid(row=3, column=0, columnspan=2, pady=(20, 10), sticky="w")

        constraint_entries = []
        senses = []
        rhs_entries = []

        for i in range(num_constraints):
            frame_label = ttk.LabelFrame(input_frame, text=f"Constraint {i+1}", padding=10)
            frame_label.grid(row=4+i, column=0, columnspan=2, pady=10, sticky="w")

            row_entries = []
            col_index = 0
            for j in range(num_vars):
                ttk.Label(frame_label, text=f"x" + to_subscript(j + 1), font=(font_type, 13)).grid(row=0, column=col_index, padx=5, pady=5, sticky="e")
                coef_entry = ttk.Entry(frame_label, width=10)
                coef_entry.grid(row=0, column=col_index+1, padx=5, pady=5)
                row_entries.append(coef_entry)
                col_index += 2

            sense_var = tk.StringVar(value="<=")
            senses.append(sense_var)

            ttk.Label(frame_label, text="RHS:", font=(font_type, 10)).grid(row=0, column=col_index+1, padx=5, pady=5, sticky="e")
            rhs_entry = ttk.Entry(frame_label, width=10)
            rhs_entry.grid(row=0, column=col_index+2, padx=5, pady=5)
            rhs_entries.append(rhs_entry)

            constraint_entries.append(row_entries)

        solve_btn = ttk.Button(input_frame, text="Solve", 
                               command=partial(process_inputs, num_vars, num_constraints, obj_entries, constraint_entries, senses, rhs_entries))
        solve_btn.grid(row=4+num_constraints, column=0, columnspan=2, pady=10)


    except ValueError as e:
        show_error(str(e))


def process_inputs(num_vars, num_constraints, obj_entries, constraint_entries, senses, rhs_entries):
    try:
        
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

        solution, obj_value = simplex(c, A, b)
        integer_solutionbb, integer_obj_valuebb = branch_and_bound(c, A, b,constraint_type)

        gc = GomoryCut(A, b, c,constraint_type)
        integer_solutiongc, integer_obj_valuegc = gc.solve()
        
        integer_solutiongc = [round(x, 2) for x in integer_solutiongc]
        integer_solutiongc = [float(x) for x in integer_solutiongc]
        integer_obj_valuegc = round(integer_obj_valuegc, 2)

        lpsol = ""
        ilpsol = ""
        ilpsolgc = ""
        for i in range(len(solution)):
            if i !=0:
                lpsol += ", "
            lpsol += "x" + to_subscript(i + 1) +"="+str(solution[i])
        for i in range(len(integer_solutionbb)):
            if i !=0:
                ilpsol += ", "
            ilpsol += "x" + to_subscript(i + 1) +"="+str(integer_solutionbb[i])
        for i in range(len(integer_solutiongc)):
            if i !=0:
                ilpsolgc += ", "
            ilpsolgc += "x" + to_subscript(i + 1)  +"="+str(integer_solutiongc[i])

        lp_solution = f"Simplex Solution: {lpsol}, Optimal Value: {obj_value}"
        ilp_solutionbb = f"Branch and Bound Solution: {ilpsol}, Optimal Value: {integer_obj_valuebb}"
        ilp_solutiongc = f"Gomory Cut Solution: {ilpsolgc}, Optimal Value: {integer_obj_valuegc}"

        show_results(lp_solution, ilp_solutionbb, ilp_solutiongc)

    except ValueError as e:
        show_error(f"Invalid input: {e}")
    except Exception as e:
        show_error(f"Error: {e}")


root = tk.Tk()


root_color = "#f2f3f5"

root.geometry("1000x800")
root.configure(bg=root_color)
root.title("Integer Linear Programming Solver")


canvas = tk.Canvas(root, bg=root_color, highlightthickness=0) 
scrollbar = ttk.Scrollbar(root, orient="vertical", command=canvas.yview)
canvas.configure(yscrollcommand=scrollbar.set)

#add another scrollbar for x-axis
scrollbarx = ttk.Scrollbar(root, orient="horizontal", command=canvas.xview)
canvas.configure(xscrollcommand=scrollbarx.set)



main_frame = tk.Frame(canvas, bg=root_color) 
canvas.create_window((0, 0), window=main_frame, anchor="nw")

scrollbar.pack(side="right", fill="y")
scrollbarx.pack(side="bottom", fill="x")
canvas.pack(side="left", fill="both", expand=True)

main_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))

# Top frame for inputs

top_frame = tk.Frame(main_frame, bg=root_color)  
top_frame.pack(side="top", fill="x", pady=10)  


top_frame.grid_columnconfigure(0, weight=1)  
top_frame.grid_columnconfigure(6, weight=1)  


tk.Label(top_frame, text="Number of Variables:", bg=root_color, font=(font_type, 14, 'bold')).grid(row=0, column=1, padx=5, pady=5, sticky="e")
entry_vars = tk.Entry(top_frame, width=10)
entry_vars.grid(row=0, column=2, padx=5, pady=5)

tk.Label(top_frame, text="Number of Constraints:", bg=root_color, font=(font_type, 14, 'bold')).grid(row=0, column=3, padx=5, pady=5, sticky="e")
entry_constraints = tk.Entry(top_frame, width=10)
entry_constraints.grid(row=0, column=4, padx=5, pady=5)




ttk.Button(top_frame, text="Generate Input Field", command=solve_problem, style="Custom.TButton").grid(row=0, column=5, padx=10, pady=5)
ttk.Button(top_frame, text="Help", command=show_help, style="Custom.TButton").grid(row=0, column=6, padx=10, pady=5)

#exit button
ttk.Button(top_frame, text="Exit", command=root.quit, style="Custom.TButton").grid(row=0, column=7, padx=10, pady=5)


#Input frame for dynamic input 
input_frame = tk.Frame(main_frame, bg=root_color) 
input_frame.pack(side="top", fill="both", expand=True, pady=10)

#Results frame
results_frame = tk.Frame(main_frame, bg=root_color, padx=20, pady=10)
results_frame.pack(side="top", fill="x")

root.mainloop()
