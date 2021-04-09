import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


'''
Программа по проверка ДЗ № 1
'''


G = float(input('Расход, кг/с = '))
T_t = float(input('Температура торможения, К = '))
p_t = float(input('Давление торможения, Па = '))
alfa_g = float(input('угол полураскрытия диффузора (alfa/2 ), град = '))
R = float(input('Газовая постоянная, Дж / (кг*К) = '))
k = float(input('Показатель адиабаты = '))
h = float(input('Высота полёта, км = ')) * 1000


pi = np.pi

alfa = alfa_g * np.pi / 180

pi_g = [0,0,0,0,0,0]
T = [0,0,0,0,0,0]
T1 = [0,0,0,0,0,0]
ro = [0,0,0,0,0,0]
S = [0,0,0,0,0,0]
r = [0,0,0,0,0,0]
d = [0,0,0,0,0,0]
V = [0,0,0,0,0,0]
a = [0,0,0,0,0,0]
d = [0,0,0,0,0,0]
p = [0,0,0,0,0,0]
lam = [0,0,1,0,0,0]
M = [0,0,0,0,0,0]
gdf_g = [0,0,0,0,0,0]


def mm(kk,rr):
    return (2 * kk / (kk + 1)) ** 0.5 * (2 / (kk + 1)) ** (1 / (kk - 1)) / (rr) ** 0.5


def p_alt(h):
    return 101325 * (1 - ((0.0065 / 288.15) * h)) ** (9.80665 / (0.0065 * 287.053))


def ro_gdf(lam):
    return (1 - (k - 1) / (k + 1) * lam ** 2) ** (1 / (k - 1))


def pi_gdf(lam):
    return (1 - (k - 1) / (k + 1) * lam ** 2) ** (k / (k - 1))


def t_gdf(lam):
    return (1 - (k - 1) / (k + 1) * lam ** 2)


def g_gdf(lam):
    return ((k + 1) / 2) ** (1/ (k - 1)) * (1 - (k - 1)/(k + 1) * lam ** 2) ** (1/(k - 1)) * lam


def v_max(Temp):
    return np.sqrt(2 * k / (k - 1) * R * Temp)


def akr(Temp):
    return np.sqrt(2 * k / (k + 1) * R * Temp)


def lam_from_pi(pi_g):
    return np.sqrt((1-pi_g ** ((k - 1) / k)) * (k+1)/(k-1))


def diam(S):
    return np.sqrt(4 * S / np.pi)


def soplo_lav():
    ro_t = p_t / (R * T_t)

    #  Параметры в критическом сечении
    T[2] = T_t * t_gdf(1)
    p[2] = p_t * pi_gdf(1)
    ro[2] = ro_t * ro_gdf(1)
    V[2] = akr(T_t)
    a_kr = akr(T_t)
    #  Площадь критического сечения
    S[2] = G / (ro[2] * a_kr)
    d[2] = diam(S[2])

    #  Параметры при срезе (сечение 5)
    pi_g[5] = p_0 / p_t
    lam[5] = lam_from_pi(pi_g[5])
    V[5] = lam[5] * V[2]
    T[5] = T_t * t_gdf(lam[5])
    ro[5] = ro_t * ro_gdf(lam[5])
    a[5] = np.sqrt(k * R * T[5])
    gdf_g[5] = g_gdf(lam[5])

    #  Площадь среза
    S[5] = S[2] / g_gdf(lam[5])
    d[5] = diam(S[5])

    #  Геометрия сопла
    l2 = (d[5] - d[2]) / (2 * np.tan(alfa))  # сверхзвуковая часть
    l1 = 2 * d[2]
    l = l1 + l2
    d[0] = d[2] * 2
    d[1] = d[2] * 1.5
    d[3] = d[2] + l2 / 3 * 2 * np.tan(alfa)
    d[4] = d[2] + 2 * l2 / 3 * 2 * np.tan(alfa)

    S[0] = np.pi / 4 * d[0] ** 2
    S[1] = np.pi / 4 * d[1] ** 2
    S[3] = np.pi / 4 * d[3] ** 2
    S[4] = np.pi / 4 * d[4] ** 2

    gdf_g[0] = S[2] / S[0]
    gdf_g[1] = S[2] / S[1]
    gdf_g[3] = S[2] / S[3]
    gdf_g[4] = S[2] / S[4]

    def f_g_gdf_0(x):
        return g_gdf(x) - gdf_g[0]

    lam[0] = float(fsolve(f_g_gdf_0, 0))
    print('Приведённая скорость 1 = ',lam[0])

    def f_g_gdf_1(x):
        return g_gdf(x) - gdf_g[1]

    lam[1] = float(fsolve(f_g_gdf_1, 0))
    print('Приведённая скорость 2 = ',lam[1])

    def f_g_gdf_3(x):
        return g_gdf(x) - gdf_g[3]

    lam[3] = float(fsolve(f_g_gdf_3, 1.9))
    print('Приведённая скорость 3 = ', lam[2])
    print('Приведённая скорость 4 = ',lam[3])

    def f_g_gdf_4(x):
        return g_gdf(x) - gdf_g[4]

    lam[4] = float(fsolve(f_g_gdf_4, 1.9))
    print('Приведённая скорость 5 = ', lam[4])
    print('Приведённая скорость 6 =', lam[5])


    # Оставшиеся параметры
    for i in range(6):
        T[i] = T_t * t_gdf(lam[i])
        a[i] = np.sqrt(k * R * T[i])
        V[i] = lam[i] * a_kr
        p[i] = p_t * pi_gdf(lam[i])

    for i in range(6):
        lam[i] = float('{:.3f}'.format(lam[i]))
        p[i] = float('{:.1f}'.format(p[i]))
        V[i] = float('{:.1f}'.format(V[i]))
        T[i] = float('{:.1f}'.format(T[i]))
        a[i] = float('{:.1f}'.format(a[i]))
        r[i] = d[i] / 2 * 1000
        M[i] = V[i] / a[i]
        M[i] = float('{:.2f}'.format(M[i]))
        r[i] = float('{:.2f}'.format(r[i]))
        d[i] = float('{:.4f}'.format(d[i]))

    print(f'a кр = {a_kr}  м/с')
    print('Ламбда = ', lam)
    print('p =  ', p, 'Па')
    print('a =  ', a, 'м/с')
    print('V =  ', V, 'м/с')
    print('М =  ', M, '')
    print('T =  ', T, 'K')
    print('r = ', r, 'мм')
    print('d = ', d, 'м')

    print('_____Длина сопла_______')
    print('Длина дозвуковой части = ', l1)
    print('Длина сверхзвуковой части =', l2)
    print('Общая длина = ', l1 + l2)

    #  Графики
    x = [1, 2, 3, 4, 5, 6]

    plt.grid()  # включение отображение сетки
    plt.xlabel("сечения, номера точек")  # ось абсцисс
    plt.ylabel("м/с")  # ось ординат

    plt.vlines(3, V[0], V[5], color='r')  # Вертикальные линии

    plt.plot(x, a, 'r--', label='скорость звука')
    plt.plot(x, V, 'g', label='скорость V')
    plt.legend()

    plt.show()

    return p, lam, S[5]


p_0 = p_alt(h)
pp, lamb, F5 = soplo_lav()
p_2 = 101325

def shock_wave(p_2=101325):

    print(f'Давление за соплом на уровне земли = {p_2} Па')

    #  Пусть скачок в 6-м сечении (на среде сопла)

    lam_za_sk = 1 / lamb[5]  # Приведённая скорость за скачок
    print(f'Приведённая скорость за скачком, если он в 6 сечении {lam_za_sk}')
    print(g_gdf(lamb[5]))
    p_t_za_sk = p_t * g_gdf(lamb[5]) / g_gdf(lam_za_sk)
    p_za_sk_6 = p_t_za_sk * pi_gdf(lam_za_sk)
    print('прив ск за скачком = ', lam_za_sk)
    print('давл торм за СК =', p_t_za_sk)
    print('давление за скачком = ', p_za_sk_6)
    if p_za_sk_6 >= p_2:
        print('Давление за скачком выше, скачок уплотнения будет на срезе сопла')
    else:
        print('Скачок в сопле')
        lam_var1, p_6_posle_ck1 = find_shock()
        return lam_var1, p_6_posle_ck1




def find_shock():
    print('Поиск скачка уплотнения...')
    lam_var = np.arange(1, lamb[5], 0.01)

    lam_za_sk_var = 1 / lam_var
    p_za_sk_var = p_t * g_gdf(lam_var) / g_gdf(lam_za_sk_var)  #  давление за скачком уплотнения в сопле

    f_var = G * np.sqrt(T_t) / (mm(k,R) * p_t * g_gdf(lam_var))
    d_var = np.sqrt(4 * f_var / pi)

    g_gdf_var_za_ck = f_var * g_gdf(lam_za_sk_var) / F5

    lam_6_posle_ck = []

    ### определение приведённой скорости через ГДФ Расхода
    for i in range(len(g_gdf_var_za_ck)):
        gdf_number = g_gdf_var_za_ck[i]
        def ff_gg(x):
            return g_gdf(x) - gdf_number
        lam_6_posle_ck.append(float(fsolve(ff_gg, 0)))



    lam_6_posle_ck = np.array(lam_6_posle_ck)

    p_6_posle_ck = p_za_sk_var * pi_gdf(lam_6_posle_ck)


    plt.grid()  # включение отображение сетки
    plt.xlabel("Приведённая скорость")  # ось абсцисс
    plt.ylabel("Давление, Па ")  # ось ординат

    #plt.vlines(3, V[0], V[5], color='r')  # Вертикальные линии
    plt.axhline(y=101325, xmin=0, xmax=2, label='Атмосфера')
    plt.plot(lam_var, p_6_posle_ck, 'r--', label='Давление на срезе сопла')
    #plt.plot(x, V, 'g', label='скорость V')
    plt.legend()

    plt.show()


    return lam_var, p_6_posle_ck

def find_intersection(t, curve1, curve2):
    '''
    Функция нахождения пересечения графиков
    :return: точка пересечения двух графиков
    '''
    intersections = []
    prev_dif = 0
    t0, prev_1, prev_2 = None, None, None

    for t1, c1, c2 in zip(t, curve1, curve2):
        new_dif = c2 - c1
        print(c1, c2)
        if np.abs(new_dif) < 0.00001:
            intersections.append((t1, c1))
        elif new_dif * prev_dif < 0:
            denom = prev_dif - new_dif
            intersections.append(((-new_dif * t0 + prev_dif * t1) / denom, (c1 * prev_c2 - c2 * prev_c1) / denom))
            print('оппа!')
            print('new_dif =', new_dif)
            print('prev_dif =', prev_dif)
            print('denom = ', denom)
            print(t0, t1)
            print(c1, c2)

        t0, prev_c1, prev_c2, prev_dif = t1, c1, c2, new_dif

    #print('Вроде нашли, но это не точно!!!')
    #print(intersections)

    return intersections

def find_geom_shock(L):
    '''
    Нахождение координаты X скачка и Диаметра, где происходит скачок
    :return: x
    '''

    x = np.arange(0, L, 0.02)

    pass



def draw_geometry():
    '''
    Рисуем геометрию спрофилированного сопла
    :return: рисунок
    '''
    

    pass


lam_var_2, p_6_posle_ck_2 = shock_wave()
p_2_massive = []
for i in range(len(lam_var_2)):
    p_2_massive.append(p_2)


print('Последние значения ')
print(lam_var_2)
print(p_6_posle_ck_2)

intersection_find = find_intersection(lam_var_2, p_6_posle_ck_2, p_2_massive)
print(intersection_find)


