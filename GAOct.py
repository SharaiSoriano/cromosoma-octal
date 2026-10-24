import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import random
import copy
from typing import List, Dict

from CromosomaOct import CromosomaOct
from interfaz import Aptitud8Reinas, SeleccionTorneo, ReproduccionUnCruce, MutacionFlip

class GA(object):
    def __init__(self, num_inv: int, num_genes: int, hyperparams: Dict):
        self.poblacion: List[CromosomaOct] = [CromosomaOct(num_genes) for i in range(num_inv)]
        self.num_inv: int = num_inv
        self.hyperparams: Dict = hyperparams
        self.evaluar = hyperparams["Aptitud"]
        self.seleccion = hyperparams["Seleccion"]
        self.reproducion = hyperparams["Reproduccion"]
        self.mutacion = hyperparams["Mutacion"]
        self.num_remplazo: int = int(hyperparams["por_remplazo"] * num_inv)

        for i in range(num_inv):
            self.poblacion[i].llenar_random_data()

    def buscar_solucion(self, tolerancia: float, var: Dict = None):
        avg_tol: float = 99999999.99
        sum_score_prev: float = 0.0
        generacion: int = 1

        fail: bool = False
        fail_count: int = 0
        fail_score: float = 0.0

        while not fail:
            self.evaluar.evaluar_poblacion(self.poblacion, var)
            sum_scores: float = sum(cromosoma.score for cromosoma in self.poblacion)
            sum_scores /= self.num_inv

            if sum_scores != 0:
                avg_tol = abs(sum_score_prev - sum_scores) / sum_scores

            if (generacion + 1) % 100 == 0:
                print(f"[{generacion:05d}] avg_err : {avg_tol:.5f} avg_score: {sum_scores:.5f}")

            mejor_score = max(c.score for c in self.poblacion)
            if mejor_score == 28.0:
                print(f"Solución óptima (0 conflictos) encontrada en la generación {generacion}")
                break

            if fail_score == avg_tol:
                fail_count += 1

            if fail_count >= 10:
                fail = True
                break

            if avg_tol < tolerancia:
                break

            fail_score = avg_tol
            sum_score_prev = sum_scores

            self.poblacion = sorted(self.poblacion, key=lambda c: c.score, reverse=True)
            mating_pool = self.seleccion.seleccion_poblacion(self.poblacion)
            generacion_nueva = copy.deepcopy(self.poblacion[: self.num_remplazo])

            self.reproducion.reproduccion_poblacion(generacion_nueva, mating_pool, self.num_inv)
            self.mutacion.mutacion_poblacion(generacion_nueva)

            generacion += 1
            self.poblacion = copy.deepcopy(generacion_nueva)

def imprimir_tablero(cromosoma: CromosomaOct):
    print("\n Tablero")
    for fila in range(8):
        linea = ""
        for col in range(1, 9):
            if cromosoma.get_bits(col) == fila:
                linea += " R "
            else:
                linea += " . "
        print(linea)


if __name__ == "__main__":
    num_genes = 8   # 8 columnas = 8 genes octales
    cant_ind = 100  # Tamaño de población

    gaparams = {}
    gaparams["Aptitud"] = Aptitud8Reinas()      
    gaparams["Seleccion"] = SeleccionTorneo(3)   
    gaparams["Reproduccion"] = ReproduccionUnCruce(num_genes)
    gaparams["Mutacion"] = MutacionFlip(0.05)
    gaparams["por_remplazo"] = 0.5

    algoritmoGA = GA(cant_ind, num_genes, gaparams)
    algoritmoGA.buscar_solucion(1e-4)

    solucion = sorted(algoritmoGA.poblacion, key=lambda c: c.score, reverse=True)
    
    print(f"\n{'id':<7} {'Cromosoma Octal':<20} {'Fitness':<10} {'Conflictos':<10}")
    for i, cromosoma in enumerate(solucion[:15]):
        conflictos = 28 - int(cromosoma.score)
        print(f"[{i + 1:>2} ] {cromosoma.a_octal_string():<20} {cromosoma.score:<10.0f} {conflictos:<10}")

    mejor = solucion[0]
    print(f"\nMejor Solución: {mejor.a_octal_string()} | Fitness: {mejor.score:.0f} | Conflictos: {28 - int(mejor.score)}")
    imprimir_tablero(mejor)