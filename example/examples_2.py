#!

import phdtools.descriptors.QUIP as QUIP

desc = QUIP.Descriptors(desc_type='soap', 
                        params_file='quip_soap.json')
print(desc.params_dict)