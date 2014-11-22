#    Author: Philipp Eller
import sys
class progbar:
    def __init__(self,width):
        self.width=width
        sys.stdout.write("[%s]" % ("-" * width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (width+1)) # return to start of line, after '['
        self.pos = 0
        self.done = False
    def move(self):
        if not self.done:
            sys.stdout.write("*")
            self.pos+=1
            sys.stdout.flush()
            if self.pos == self.width:
                sys.stdout.write("]\n")
                sys.stdout.flush()
                self.done = True


if __name__ == "__main__":
    from time import sleep
    pb = progbar(10)
    for i in range(10):
        sleep(0.1)
        pb.move()
