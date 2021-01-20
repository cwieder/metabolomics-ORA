import xml.etree.ElementTree as ET

root = ET.parse('/Users/cw2019/Downloads/csf_metabolites.xml').getroot()

for type_tag in root.findall('kind'):
    value = type_tag.get('kind')
    print(value)