from timeit import timeit

import sys

def timed(func):
    def wrapped(*args, **kwargs):
        output = [None]
        def saveoutput():
            output[0] = func(*args, **kwargs)
        time = timeit(saveoutput, number=1)
        return output[0], time
    return wrapped

def main():
    @timed
    def timed_tester(n):
        return len(range(n))

    print timed_tester(5)
    print timed_tester(10000000)

    return 0

if __name__ == '__main__':
    sys.exit(main())
