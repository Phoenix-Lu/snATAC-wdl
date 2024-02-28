import os
import json

def Save_file_Path(json_path, path_key, path_dicts:dict):
    with open(json_path, "r") as json_file:
        data = json.load(json_file)
        for file_name, file_path in path_dicts.items():
            if path_key not in data:
                data[path_key] = {}
            
            data[path_key][file_name] = file_path

    with open(json_path, "w") as json_file:
        json.dump(data, json_file)
    return json_path






def main():
    Save_file_Path('/share/home/lushaorong/test/test.json',
    'class3', {"file3":"file3path"})
    
if __name__ == '__main__':
    main()



    
