from cx_Freeze import setup, Executable
import matplotlib
import reactor
import os.path
PYTHON_INSTALL_DIR = os.path.dirname(os.path.dirname(os.__file__))
os.environ['TCL_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tcl8.6')
os.environ['TK_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tk8.6')

build_exe_options = {"includes":["matplotlib.backends.backend_tkagg"],
                     "include_files":[(matplotlib.get_data_path(), "mpl-data"),os.path.join(PYTHON_INSTALL_DIR, 'DLLs', 'tk86t.dll'),
            os.path.join(PYTHON_INSTALL_DIR, 'DLLs', 'tcl86t.dll')],
                     "excludes":[], "optimize":2, "includes":['reactor'],
                     }

setup(name = "Reactor simulation" ,
      version = "0.1" ,
      options = {"build_exe": build_exe_options},
      description = "Reactor simulation" ,
      executables = [Executable("NuclearReactorSimulator.py")])