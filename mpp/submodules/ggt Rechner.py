# -*- coding: utf-8 -*-
"""
Created on Sat Jun 29 10:47:47 2024

@author: phili
"""

import random

def ggt(a, b):
    m = max(a, b)
    n = min(a, b)
    
    while n != 0:
        m, n = n, m % n
    
    return m



def eea(a, b):
    if b == 0:
        return (1, 0)
    
    x1, x2, y1, y2 = 1, 0, 0, 1
    
    while b != 0:
        q = a // b
        r = a % b
        a, b = b, r
        x = x1 - q * x2
        y = y1 - q * y2
        x1, x2 = x2, x
        y1, y2 = y2, y
    
    return (x1, y1)


def test_random_pairs():
    for _ in range(100):
        a = random.randint(1, 100)
        b = random.randint(1, 100)
        
        g = ggt(a, b)
        x, y = eea(a, b)
        
        print(f"Für a = {a}, b = {b}:")
        print(f"  ggT({a}, {b}) = {g}")
        print(f"  Lösung der Gleichung {a}*{x} + {b}*{y} = {g}")
        print(f"  Überprüfung: {a}*{x} + {b}*{y} - {g} = {a*x + b*y - g}")
        print()  # Leerzeile zur besseren Lesbarkeit

test_random_pairs()