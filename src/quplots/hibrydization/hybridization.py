from .sp import sp
from .sp2 import sp2
from .sp3 import sp3
from .sp2d import sp2d
from .sp3d import sp3d
from .sp3d2 import sp3d2
from .sp3d import sp3d  

class hybridization:
    def __init__(self):
        pass

    def generate_sp(self):
        return sp()

    def generate_sp2(self):
        return sp2()

    def generate_sp3(self):
        return sp3()

    def generate_sp2d(self):
        return sp2d()

    def generate_sp3d(self):
        return sp3d()

    def generate_sp3d2(self):
        return sp3d2()

    def generate_sp3d(self):
        return sp3d()
