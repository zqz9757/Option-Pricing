
#二分法 bisection method
#法一：
def bisection(func, a, b, tol=0.1, maxiter=10):
    c = (a+b)*0.5  # Declare c as the midpoint ab
    n = 1  # Start with 1 iteration
    while n <= maxiter:
        c = (a+b)*0.5
        if func(c) == 0 or abs(a-b)*0.5 < tol:
            # Root is found or is very close
            return c, n
        n += 1
        if func(c) < 0:
            a = c
        else:
            b = c
    return c, n
root, iterations = bisection(lambda x: diff_K(x)-y, -5, 5, 0.00001, 100)
print("Root is:",root)
print("Iterations:",iterations)
#法二:
def bisection_cal(func,start,end,precision):
    d=func(start)
    u=func(end)
    if d*u>0:
        print('No solution')
        return False
    if d==0:
        print('The solution is:',start)
        return start
    if u==0:
        print('The solution is:',end)
        return end
    while abs(start-end)>precision:
        mid=(start+end)/2
        m=func(mid)
        if m==0:
            print('The solution is:',mid)
            return mid
        if d*m<0:
            end=mid
        if m*u<0:
            start=mid
    print('The solution is:',end)
    return end

