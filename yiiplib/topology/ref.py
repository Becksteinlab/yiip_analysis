from pkg_resources import resource_filename
from .universe import YiiPUniverse
import os

ref_if = resource_filename(__name__, '5VRF_holo.pdb')
ref_of = resource_filename(__name__, '5VRF_of.pdb')
ref_if_apo = resource_filename(__name__, '7KZX_apo.pdb')
ref_if_holo = resource_filename(__name__, '7KZZ_holo.pdb')

ref_if = YiiPUniverse(ref_if)
ref_of = YiiPUniverse(ref_of)
ref_if_apo = YiiPUniverse(ref_if_apo)
ref_if_holo = YiiPUniverse(ref_if_holo)