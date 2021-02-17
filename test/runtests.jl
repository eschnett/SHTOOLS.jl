using SHTOOLS
using Test

p = PlmBar(4, 0)
p, dp = PlmBar_d1(4, 0)
p = PlBar(4, 0)
p, dp = PlBar_d1(4, 0)

p = PlmON(4, 0)
p, dp = PlmON_d1(4, 0)
p = PlON(4, 0)
p, dp = PlON_d1(4, 0)

p = PlmSchmidt(4, 0)
p, dp = PlmSchmidt_d1(4, 0)
p = PlSchmidt(4, 0)
p, dp = PlSchmidt_d1(4, 0)

p = PLegendreA(4, 0)
p, dp = PLegendreA_d1(4, 0)
p = PLegendre(4, 0)
p, dp = PLegendre_d1(4, 0)

index = PlmIndex(4, 3)
@test index == 14
