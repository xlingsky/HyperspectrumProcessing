import os
import sys
import json
import xml.etree.ElementTree as ET

def xml_to_dict(element : ET):
    result = {}
    
    if element.attrib:
        result['<xmlattr>'] = element.attrib
        
    if element.text and element.text.strip():
        result = element.text.strip()
        
    for child in element:
        child_data = xml_to_dict(child)

        if child.tag in result:
            if isinstance(result[child.tag], list):
                result[child.tag].append(child_data)
            else:
                result[child.tag] = [result[child.tag], child_data]
        else:
            result[child.tag] = child_data

    return result

def xml_string_to_json(xml_string):
    root = ET.fromstring(xml_string)
    data_dict = xml_to_dict(root)
    return json.dumps({root.tag:data_dict}, indent=2, ensure_ascii=False)

# def generate_object_detection_xml(json : dict):

# def json_to_tasks(json : str, taskdir : os.path):

# def object_detection():

if __name__ == "__main__":
    import sys
    with open(sys.argv[1], 'r') as f:
        xml_string = f.read()
        with open(sys.argv[2], 'w') as f:
            f.write(xml_string_to_json(xml_string))