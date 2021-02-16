#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@文件        :module_completeness.py
@说明        :
@时间        :2020/12/30 15:27:14
@作者        :Lv Yongxin
@版本        :1.0
'''


class Element():
    def __init__(self, express="", additional_info=""):
        # print(express)
        self.is_chain = True
        self.elements = []
        self.ko = ""
        self.additional_info = additional_info
        self.__calculate(express)

    def __calculate(self, express: str):
        """Calculate express"""
        if not express:
            return
        p = len(express) -1
        if p == 5:  # ko
            self.ko = express
            return
        no_comma = []
        while p >= 0:
            c = express[p]
            if c in "0123456789":  # get a brief KO numbner
                no_comma.append(Element(express[p - 5:p + 1]))
                p -= 5
            elif c in ")":  # you are trapped into a sub element
                last_p = p
                bracket_stack = 1
                plb = express.rfind("(", 0, p)
                prb = express.rfind(")", 0, p)
                while bracket_stack:
                    if plb < prb:
                        bracket_stack += 1
                        p = prb
                        prb = express.rfind(")", 0, p)
                    else:
                        bracket_stack -= 1
                        p = plb
                        plb = express.rfind("(", 0, p)
                no_comma.append(Element(express[p + 1:last_p]))
            elif c in ",":
                if self.is_chain:
                    self.is_chain = False
                self.elements.append(Element.elist(no_comma))
                no_comma = []
                # Stupid KEGG think "(K09880,K08965 K08966)" is equal to "(K09880,(K08965 K08966))".
            p -= 1
        if self.elements:
            self.elements.append(Element.elist(no_comma))
        else:
            self.elements = no_comma
        # print(*self.elements)
        self.elements.reverse()

    def completeness(self, ko_match: dict):
        """Complessness of given match, ko is its dict"""
        count = 0
        if self.ko:
            return 1 if self.ko in ko_match else 0
        # multipy elements
        if self.is_chain:
            for element in self.elements:
                count += element.completeness(ko_match)
            return count / len(self.elements)
        # self.is_chain is False
        return max(
            [element.completeness(ko_match) for element in self.elements])

    def __str__(self):
        if self.ko:
            return self.ko
        sep = " " if self.is_chain else ","
        return sep.join([
            kid.ko if kid.ko else "(" + str(kid) + ")" for kid in self.elements
        ])

    @classmethod
    def elist(cls, no_comma, is_chain=True, additional_info=""):
        """New Element by a list"""
        if len(no_comma) == 1:
            e = no_comma[0]
        else:
            e = Element()
            e.elements = no_comma
            e.is_chain = is_chain
            e.additional_info = additional_info
        return e


dict_module = {}
#fin = open("/lustre/home/acct-ioozy/ioozy-user2/LYX/Script/pipe_temp/module_KO_KEGG.txt").read().split('\n')[:-1]
fin=open("def_module_ko.txt").read().split('\n')[:-1]
ptr = 0
#title = "LYX"
while ptr < len(fin):
    if fin[ptr].strip()[0:5] == "ENTRY":
        module_name = fin[ptr].split()[1]
#        print(fin[ptr])
        ptr += 1
        line = " ".join(fin[ptr].strip().split()[1:])
        name = line.strip()
        ptr += 1
        module = " ".join(fin[ptr].strip().split()[1:])
        dict_module[module_name+'\t'+name] = Element(module, module_name)
#    else:  # title
#        title = fin[ptr].strip()[:-1]
#        if title:
#            dict_module[title] = {}
    ptr += 1
#print(dict_module)

def kegg_module(ko_match):
    """As you like"""
    out_data = {}
    metabolism_data = {}
    for metabolism, element in dict_module.items():
        value = element.completeness(ko_match)
        if value:
#                print(element.additional_info)
                metabolism_data[metabolism] = value
        out_data.update(metabolism_data)
    return out_data

"""
    for pathway in dict_module:
        metabolism_data = {}
        for metabolism, element in dict_module[pathway].items():
            value = element.completeness(ko_match)
            if value:
                print(element.additional_info)
                metabolism_data[metabolism] = value
        out_data.update(metabolism_data)
    return out_data
"""

from sys import argv
sc, infile, outfile = argv
dict_input = {}
for line in open(infile):
    temp = line.strip().split('\t')
    gene = temp[0]
    ko = temp[1]
    dict_input[ko] = 1
#print(dict_input)
#print(kegg_module({"K22480":1}))
#print(kegg_module({"K14080":1,"K00125":1,"K22480":1}))

dict_out = kegg_module(dict_input)
#print(dict_out)
fout = open(outfile, 'w')
for key in dict_out.keys():
    fout.write(key + '\t' + str(dict_out[key]) + '\n')

