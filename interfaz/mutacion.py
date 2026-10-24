import random
from abc import ABC, abstractmethod
from typing import List
from CromosomaOct import CromosomaOct

class Mutacion(ABC):
    @abstractmethod
    def mutacion_poblacion(self, poblacion: List[CromosomaOct]):
        pass

class MutacionFlip(Mutacion):
    def __init__(self, probabilidad: float):
        super().__init__()
        self.prob = probabilidad

    def mutacion_poblacion(self, poblacion: List[CromosomaOct]):
        max_valor = 2**poblacion[0].gen_width - 1 # 7
        for crom in poblacion:
            for i in range(1, crom.n_genes + 1):
                prob_gen = random.random()
                if prob_gen < self.prob:
                    
                    nuevo_valor = random.randint(0, max_valor)
                    crom.set_bits(nuevo_valor, i)