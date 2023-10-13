class FM:
    def __init__(self):
        self.model = True


class RR:
    def __init__(self):
        self.model = "a"


class MODELS:
    def __init__(self) -> None:
        self.fm = FM()
        self.rr = RR()
        print(self.fm.model)
        self.rr.fm = self.fm
        self.rr.fm.model = False
        print(self.fm.model)
        print(self.rr.fm.model)
        self.fm.model = True
        print(self.fm.model)
        print(self.rr.fm.model)


models = MODELS()
