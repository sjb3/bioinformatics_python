def fib(num):
    old, new = 1, 1
    for i in range(num-1):
        temp = new
        new = old
        old = temp + new

    return new


print(fib(12))


def fib_ht(num, ht={1: 1, 2: 1, 3: 2}):
    if num in ht:
        return ht[num]

    ht[num] = fib_ht(num-1, ht) + fib_ht(num-2, ht)
    return ht[num]


print(fib_ht(12))


def rabbit_recur(mon, off):
    if mon < 2:
        return mon
    else:
        return rabbit_recur(mon-1, off) + rabbit_recur(mon-2, off) * off


print(rabbit_recur(5, 2))
'''
5 months x 2 fecundity

Month 1: [o]
Month 2: [O]
Month 3: [O] [o] [o]
Month 4: [O] [o] [o] [O] [O]
Month 5: [O] [o] [o] [O] [O] [o] [o] [o] [o]
Month 6: [O] [O] [O] [o] [o] [o] [o] [o] [o] [o] [o]
'''
