def kalman_factory(class_name, options):
    import importlib
    import os
    import sys

    # Ensure the current directory is included in sys.path for module lookup
    current_directory = os.path.dirname(__file__)
    if current_directory not in sys.path:
        sys.path.append(current_directory)

    # Split the class name to determine the module it belongs to
    if 'KalmanConventional' in class_name:
        module_name = 'kalman_conventional'  # File name without .py
    elif 'KalmanUltimate' in class_name:
        module_name = 'kalman_ultimate'  # File name without .py
    else:
        raise ValueError(f"Module for class {class_name} not found.")

    # Import the module dynamically
    try:
        module = importlib.import_module(module_name)  # Import the correct module
        class_ = getattr(module, class_name)  # Get the class by name
        return lambda: class_(options)  # Return a callable to create the class instance
    except AttributeError:
        raise ValueError(f"Class {class_name} not found in the {module_name} module.")
