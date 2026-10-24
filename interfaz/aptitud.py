from abc import ABC, abstractmethod
from typing import List, Dict
from CromosomaOct import CromosomaOct

class Aptitud(ABC):
    @abstractmethod
    def evaluar_poblacion(self, poblacion: List[CromosomaOct], var: Dict = None):
        """interfaz encargada para darle un score a la población"""
        pass

class Aptitud8Reinas(Aptitud):
    def evaluar_poblacion(self, poblacion: List[CromosomaOct], var: Dict = None):
        for cromosoma in poblacion:
            conflictos: int = 0
            n: int = cromosoma.n_genes

            for i in range(1, n + 1):
                fila1 = cromosoma.get_bits(i)
                for j in range(i + 1, n + 1):
                    fila2 = cromosoma.get_bits(j)
                    
                    if fila1 == fila2:
                        conflictos += 1
                    elif abs(fila1 - fila2) == abs(i - j):
                        conflictos += 1

            cromosoma.score = float(28 - conflictos)