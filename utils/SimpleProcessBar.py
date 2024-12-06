class SimpleProcessBar:
    def __init__(self, total, finished='■', unfinished='□'):
        self.total = total
        if self.total > 100:
            self.ratio = 100 / self.total 
            self.total = 100
            self.finished = 1
        else:
            self.ratio = 1
            self.finished = 0
        # self.finished = 0
        self.unfinished_char = unfinished
        self.finished_char = finished
        self.bar = [self.unfinished_char] * 100

    def generate_bar(self):
        finished = int(self.finished / self.total * 100)
        self.bar[:finished] = self.finished_char * finished
        return ''.join(self.bar) + f' {finished} %'

    def incr(self):
        self.finished += self.ratio
        if self.finished > 100:
            self.finished = 100
        bar = self.generate_bar()
        print(bar, end='\r')
