import random
from typing import Dict, List

class CromosomaOct(object):
    def __init__(self, genes: int = 8):
        self.n_genes: int = genes
        self.score: float = -1.0
        
        self.gen_width: int = 3
        self.mascaras: Dict = {}
        self.bits: int = 0x00000000

        
        max_value: int = 2**self.gen_width - 1 # 7
        for g, i in enumerate(reversed(range(0, 24, self.gen_width))):
            self.mascaras["p" + str(g + 1)] = i
            self.mascaras["m" + str(g + 1)] = max_value << i

    def llenar_random_data(self):
        max_value = 2**self.gen_width - 1 
        for i in range(self.n_genes):
            valor = random.randint(0, max_value)
            self.set_bits(valor, i + 1)

    def registrar_mascaras(self, mascara: int, gen: int) -> None:
        self.mascaras["m" + str(gen)] = mascara

    def get_bits(self, gen: int) -> int:
        return (self.bits & self.mascaras["m" + str(gen)]) >> self.mascaras["p" + str(gen)]

    def set_bits(self, valor: int, gen: int) -> None:
        self.bits = (self.bits & ~self.mascaras["m" + str(gen)]) | (valor << self.mascaras["p" + str(gen)])

    def a_octal_string(self) -> str:
        return "".join(str(self.get_bits(i + 1)) for i in range(self.n_genes))