
a = 10.0
b = 14.87

a_0 = a*b+a*s*b+2*s*b
a_1 = a*s + 2*s*b
a_2 = a+a*s+2*s+1
a_3 = a_2*a_1 - a_0

print("a_0: {a_0}")
print("a_1: {a_1}")
print("a_2: {a_2}")
print("a_3: {a_3}")
