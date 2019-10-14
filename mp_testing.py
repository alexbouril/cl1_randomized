# # Python program to illustrate
# # the concept of locks
# # in multiprocessing
# import multiprocessing
#
#
# # function to withdraw from account
# def withdraw(balance, lock):
#     for _ in range(10000):
#         lock.acquire()
#         balance.value = balance.value - 1
#         lock.release()
#
#
# # function to deposit to account
# def deposit(balance, lock):
#     for _ in range(10000):
#         lock.acquire()
#         balance.value = balance.value + 1
#         lock.release()
#
#
# def perform_transactions():
# 	# initial balance (in shared memory)
#     balance = multiprocessing.Value('i', 100)
# 	# creating a lock object
#     lock = multiprocessing.Lock()
# 	# creating new processes
#     p1 = multiprocessing.Process(target=withdraw, args=(balance,lock))
#     p2 = multiprocessing.Process(target=deposit, args=(balance,lock))
# 	# starting processes
#     p1.start()
#     p2.start()
#
# 	# wait until processes are finished
#     p1.join()
#     p2.join()
#
# 	# print final balance
#     print("Final balance = {}".format(balance.value))
#
# if __name__ == "__main__":
#     for _ in range(10):
# 		# perform same transaction process 10 times
#         perform_transactions()

from functools import lru_cache
import multiprocessing
import os
import numpy as np


def square(n):
    # print("Worker process id for {0}: {1}".format(n, os.getpid()))
    s = 0
    for i in range(1000):
        s+=i
    return (n * n)


def inner(n):
    s=0
    for i in range(1000):
        s+=i
    return n*n


def squares_1_to_(n):
    li = np.arange(10000)
    p = multiprocessing.Pool()
    result = p.map(inner, li)
    return result



if __name__ == "__main__":
    # input list
    import datetime
    import timeit
    a = datetime.datetime.now()
    mylist = np.arange(10000)
    # creating a pool object
    p = multiprocessing.Pool()
    # map list to target function
    result = p.map(square, mylist)
    b = datetime.datetime.now()
    print(b-a)
    a = datetime.datetime.now()
    result = [0]*10000
    for i in range(10000):
        result[i]=square(i)
    b = datetime.datetime.now()
    print(b-a)
    a = datetime.datetime.now()

    r = squares_1_to_(10000)
    b = datetime.datetime.now()
    print(b-a)
