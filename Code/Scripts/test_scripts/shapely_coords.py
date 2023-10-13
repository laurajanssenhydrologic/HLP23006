from shapely.geometry import LineString, Point

a = Point(1, 2)
b = Point(2, 3)
ab = LineString([a, b])
print(ab.coords[:])

for x, y in ab.coords[:]:
    print(x, y)
