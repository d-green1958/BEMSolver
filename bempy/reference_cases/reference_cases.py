def ref_case_path(ref_name):
    import os
    from bempy.exceptions import InputDataError
    file_name = ref_name+".case"
    file_path = os.path.join(os.path.dirname(__file__), file_name)
    if os.path.exists(file_path):
        return file_path
    else:
        print("PATH:", file_path)
        raise InputDataError("Reference case not found! For a list of reference cases use bempy.list_ref_cases().")
    
    
def list_ref_cases():
    import os
    print("Reference cases are:")
    
    files = os.listdir(os.path.dirname(__file__))
    for file in files:
        if file.endswith(".case"):
            print(f"-- {file.split('.case')[0]}")
    
    
    