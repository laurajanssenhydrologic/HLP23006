class A:
    pass


class B:
    attribute = True


a = A()
a.b = B()

print(a.b.attribute)

string = "b.attribute"
ss = string.split(".")

obj = a
for s in ss:
    obj = getattr(obj, s)
    print(obj)
