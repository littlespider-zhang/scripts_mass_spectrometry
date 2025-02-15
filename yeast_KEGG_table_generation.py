import json
import requests

def up2kegg(og='sce'):

    # UniProt ID  KEGG Gene ID
    print(f"Converting UniProt ID to KEGG ID for organism {og} ...", end='')
    url = f"https://rest.kegg.jp/conv/{og}/uniprot"
    lines = requests.get(url).text.split('\n')
    # Generate the id map
    map = dict()
    for l in lines:
        if l != '':
            items = l.split('\t')
            up = items[0][3:].strip()
            name = items[1].strip()

            map[up] = name

    print(f"\tSuccess!")

    test = list(map.keys())[0]
    print(f"\teg: {test} -> {map[test]}")

    return map

def kegg2module(og='sce'):

    # KEGG ID  Module ID
    print(f"Converting KEGG ID to Module ID for organism {og} ...", end='')
    url = f"https://rest.kegg.jp/link/module/{og}"
    lines = requests.get(url).text.split('\n')

    # Generate the id map
    map = dict()
    for l in lines:
        if l != '':
            items = l.split('\t')
            up = items[0].strip()
            name = items[1][7:].strip()

            # Multiple assignments occurs
            try:
                check = map[up]
                map[up] += f",{name}"
            except:
                map[up] = name

    print(f"\tSuccess!")

    test = list(map.keys())[0]
    print(f"\teg: {test} -> {map[test]}")

    return map

def module2name():
    # KEGG ID  Module ID
    print(f"Converting Module ID to Module Name for organism ...", end='')
    url = f"https://rest.kegg.jp/list/module"
    lines = requests.get(url).text.split('\n')

    # Generate the id map
    map = dict()
    for l in lines:
        if l != '':
            items = l.split('\t')
            up = items[0].strip()
            name = items[1].strip()

            # Multiple assignments occurs
            try:
                check = map[up]
                map[up] += f",{name}"
            except:
                map[up] = name

    print(f"\tSuccess!")

    test = list(map.keys())[0]
    print(f"\teg: {test} -> {map[test]}")

    return map


def main(og='sce'):

    # Basic map: using KEGG API
    u2k = up2kegg(og=og)
    k2m = kegg2module(og=og)
    m2n = module2name()

    # Mapping UniProt ID to Module ID
    map = dict()
    for up in u2k.keys():
        map[up] = dict()
        try:
            md = k2m[u2k[up]]
            mdname = ','.join([m2n[id.strip()] for id in md.split(',')])
        except:
            print(f"Error: {up} not assigned to KEGG modules ...")
            md = ''
            mdname = ''

        map[up]["Module_ID"] = md
        map[up]["Module_Name"] = mdname

    # Writing
    with open(f"1_{og}_KEGG_module_table.json", 'w', encoding='utf-8') as f:
        json.dump(map, f, indent=2)

    return u2k,k2m,map


if __name__ == "__main__":
    # test
    u2k,k2m,map = main()
