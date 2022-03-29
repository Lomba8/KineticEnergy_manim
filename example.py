from manim import *
import math
import numpy as np
import scipy.constants as sc

c = sc.c / 10**7
gamma = MathTex(r"\gamma = 1/\sqrt{1 - v^2/c^2}")
p = MathTex(r"\mathbf {p} = m \gamma \mathbf {v}")

# L = F*dx = F*v*dt = dp/dt*v*dt = v*dp = v*d(m*v)
f1 = MathTex(
    r"d\!L = {F} \cdot d\mathbf {x} =",  # F*dx
    r"\mathbf {F} \cdot \mathbf {v} dt =",  # F*v*dt
    r"{\frac {d\mathbf {p} }{dt}}\cdot \mathbf {v} \cdot dt =",  # dp/dt*v*dt
    r"\mathbf {v} \cdot d\mathbf {p} =",  # v*dp
    r"\mathbf {v} \cdot d(m\mathbf {v} ) =",  # v*d(m*v)
    r"\mathbf {v} \cdot d(m\mathbf {v} ) =",  # L = v*d(m*v)
    substrings_to_isolate="=",
)

# L = v(dm*v + m*dv) = v*m*dv
f2 = MathTex(
    # v(dm*v + m*dv)
    r"= \mathbf {v}({d{m}\cdot \mathbf {v} } + m\cdot{d\mathbf {\mathbf {v} } })",
    r"= {m \cdot \mathbf {v} }  \cdot{d\mathbf {\mathbf {v} } }",  # v*m*dv
    substrings_to_isolate="=",
)

# d(v*v) = dv*v + v*dv = 2(v*dv)
f3 = MathTex(
    r"\text{d}( \mathbf {v} \cdot \mathbf {v} ) =",  # d(v*v)
    # dv*v + v*dv
    r"(d\mathbf {v} )\cdot \mathbf {v} +\mathbf {v} \cdot (d\mathbf {v} ) =",
    r"2(\mathbf {v} \cdot d\mathbf {v} ) ",  # 2(v*dv)
    substrings_to_isolate="=",
)

# L = m/2d(v^2) = Ek
f4 = MathTex(
    r"d\!L ="  # L=
    r"d(\frac{m}{2} {\text{d}\mathbf {v^2}}) ",  # m/2d(v^2)
    substrings_to_isolate="=",
)
# Ek = ∫F*dx = ∫v*d(m*v) = ∫d(m/2*v^2) = m*v^2/2
f4a = MathTex(
    r"E_{\kappa} =",  # Ek
    r"\int \mathbf {F} \cdot d\mathbf {x} =",  # ∫F*dx
    r"\int \mathbf {v} \cdot d(m\mathbf {v} )=",  # ∫v*d(m*v)
    r"\int d\left({\frac {mv^{2}}{2}}\right)=",  # ∫d(m/2*v^2)
    r"{\frac {m \mathbf v^{2}}{2}}",  # m*v^2/2
    substrings_to_isolate="=",
)

# L = v*d(p0) = Ek
f5 = MathTex(
    r"d\!L = \mathbf {v} \cdot d\mathbf {p} ",
    r"= E_{\text{k}",
    substrings_to_isolate="=",
)

# Ek = ∫(F*dx) = ∫v*(m*γ*v) = m*γ*v - ∫m*γ*v * dv = m*γ*v^2 - m/2*∫(γ*dv^2)
f6 = MathTex(
    r"E = \int d\!L =",  # EK = ∫L
    r"\int \mathbf {v} \cdot d\mathbf {p} =",  # ∫(v*dp)
    r"\int \mathbf {v} \cdot d(m\gamma \mathbf {v} )=",  # ∫v*(m*γ*v)
    # m*γ*v - ∫m*γ*v * dv
    r"m\gamma \mathbf {v} \cdot \mathbf {v} -\int m\gamma \mathbf {v} \cdot d\mathbf {v} =",
    # m*γ*v^2 - m/2*∫(γ*dv^2)
    r"m\gamma v^{2}-{\frac {m}{2}}\int \gamma d\left(v^{2}\right) =",
    # m*γ*v^2 - (-m*c^2/2)*∫γ*d(1 - v^2/c^2)
    r"m\gamma v^{2}-{\frac {-mc^{2}}{2}}\int \gamma d\left(1-{\frac {v^{2}}{c^{2}}}\right) =",
    # m*γ*v^2 + m*c^2*(1 - v^2/c^2)^1/2 - E0
    r"m\gamma v^{2}+mc^{2}\left(1-{\frac {v^{2}}{c^{2}}}\right)^{\frac {1}{2}}-E_{0} =",
    # m*γ*( v^2 + c^2*(1 - v^2/c^2)) - E0
    r"m\gamma \left(v^{2}+c^{2}\left(1-{\frac {v^{2}}{c^{2}}}\right)\right)-E_{0} =",
    # m*γ*( v^2 + c^2 - v^2 ) - E0
    r"m\gamma \left(v^{2}+c^{2}-v^{2}\right)-E_{0} =",
    # m*γ*c^2 - E0
    r"m\gamma c^{2}-E_{0} =",
    # (γ -1)* m*c^2
    r"(\gamma -1) mc^{2} ",
    substrings_to_isolate="=",
)

f7 = MathTex(
    r"{E_{\kappa}} = \frac{m \mathbf v^{2}}{2}",  # Ek = m/2*v^2
    r"2 \cdot {E_{\kappa}}= {m \mathbf v^{2}}",  # 2Ek=m/2*v^2
    r" \frac{2 {E_{\kappa}}}{m}= {\mathbf v^{2}}",  # 2Ek/m = v^2
    r"\mathbf {v} = \sqrt{\frac{2 \cdot {E_{\kappa}}}{m}}",  # v = √(2*ek)/m
    # lim Ek->∞ di v(Ek) =  √(2EK/m)
    r"{ \lim _{E_{\kappa}\to +\infty } \mathbf v_{E_{\kappa}} = \sqrt{\frac{2 {E_{\kappa}}}{m}}}",
    # lim Ek->∞ di v(Ek) =  √(2(+ ∞)/m)
    r"{ \lim _{E_{\kappa}\to +\infty } \mathbf v_{E_{\kappa}} = \sqrt{\frac{2 (+\infty)}{m}}}",
    # lim Ek->∞ di v(Ek) -> + ∞
    r"{ \lim _{E_{\kappa}\to +\infty } \mathbf v_{E_{\kappa}} \rightarrow +\infty}",
    # substrings_to_isolate="=",
)

f8 = MathTex(
    # Ek = (γ -1)* m*c^2
    r"{E_{\kappa}} =(\gamma -1) mc^{2}",  
    # mc^2/√(1-v^2/c^2) -mc^2
    r"{E_{\kappa}} = \frac{mc^{2}}{\sqrt{1-\frac{\mathbf v^{2}}{c^{2}}}}-mc^{2}",
    # Ek + mc^2 = mc^2*γ
    r"{E_{\kappa}} + mc^2 = mc^2\gamma",
    # Ek/mc^2 + 1 = γ
    r"\frac{E_{\kappa}}{mc^2} + 1 = \gamma",
    # Ek/mc^2 + 1 = 1/√(1-v^2/c^2)
    r"\frac{E_{\kappa}}{mc^2} + 1 = \frac{1}{\sqrt{1-\frac{\mathbf v^{2}}{c^{2}}}}",
    # [(Ek+mc^2)/mc^2]^-1 = √(1-v^2/c^2)
    r"(\frac{E_{\kappa} + mc^{2}}{mc^{2}})^{-1}= \sqrt{1-\frac{\mathbf v^{2}}{c^{2}}}",
    # mc^2/(Ek = mc^2) = √(1-v^2/c^2)
    r"\frac{mc^{2}}{E_{\kappa} + mc^{2}} = \sqrt{1-\frac{\mathbf v^{2}}{c^{2}}}",
    # - m^2*c^4/(Ek^2 + m^2*c^4) +1 = v^2/c^2
    r"- \frac{m^{2}c^{4}}{E_{\text{k}^2} + m^2c^{4}} + 1= \frac{\mathbf v^2}{c^2}",
    # [-m^2c^4 + Ek^2 + m^2c^4]/[Ek^2 +m^2c^4] = v^2/c^2
    r"\frac{-m^2c^4+E_{\kappa}^2+m^2c^4}{E_{\kappa}^2+m^2c^4} = \frac{\mathbf v^2}{c^2}",
    # Ek^2/[Ek^2 +m^2c^4] = v^2/c^2
    r"\frac{E_{\kappa}^2}{E_{\kappa}^2+m^2c^4} = \frac{\mathbf v^2}{c^2}",
    # c^2 * Ek^2/[Ek^2 +m^2c^4] = v^2
    r"c^2\frac{E_{\kappa}^2}{E_{\kappa}^2+m^2c^4} = {\mathbf v^2}",
    # c √(Ek^2/[Ek^2 +m^2c^4]) = v
    r"c\sqrt{\frac{E_{\kappa}^2}{E_{\kappa}^2+m^2c^4}} = {\mathbf v} ",
    # lim Ek->∞ di v(Ek) =  c^2 √(Ek^2/[Ek^2 +m^2c^4])
    r"{ \lim _{E_{\kappa}\to +\infty } \mathbf v_{E_{\kappa}}} = c \sqrt{\frac{E_{\kappa}^2}{E_{\kappa}^2+m^2c^4}} ",
    # c √ {Ek^2/Ek^2}* {1+[1+m^2c^4/Ek^2]}
    r" \lim _{E_{\kappa}\to +\infty } \mathbf v_{E_{\kappa}}} = c\sqrt{\frac{E_{\kappa}^2}{E_{\kappa}^2}\frac{1}{1+ \frac{m^2c^4}{E_{\kappa}^2}}}",
    # c √ {Ek^2/Ek^2}* {1+[1+m^2c^4/∞}
    r" \lim _{E_{\kappa}\to +\infty } \mathbf v_{E_{\kappa}}} = c\sqrt{\frac{1}{1+ \frac{m^2c^4}{+\infty}}}",
    # c √ {Ek^2/Ek^2}* {1+[1+0}
    r" \lim _{E_{\kappa}\to +\infty } \mathbf v_{E_{\kappa}}} = c\sqrt{\frac{1}{1+ 0}}",
    # c =  3*10^8 m/s
    r" \lim _{E_{\kappa}\to +\infty } \mathbf v_{E_{\kappa}}} = c = 3 \cdot 10^8 \frac{m}{s} ",
)

f9 = MathTex(
    # √1+x
    r"\sqrt{1+x} =",
    # (1+x)^1/2
    r"(1+x)^{\frac{1}{2}} =",
    # 1 
    r"1 ",
    # + x/2 
    r" + {\frac{x}{2}}",
    # - x^2/8
    r"- {\frac {x^{2}}{8}}",
    # + x^3/16
    r"+ {\frac {x^{3}}{16}}",
    # - 5*x^4/128
    r"- {\frac {5}{128}}x^{4}",
    substrings_to_isolate="="
)

f10 = MathTex(
    # (1+x)^alpha = 1 + alpha x + (alpha(alpha-1))/2 x^2 + (alpha(alpha-1)(alpha-2))/6 x^3 +...+ ((alpha),(n))x^n 
    r"(1+x)^\alpha = 1 + \alpha x + \frac{\alpha (\alpha -1)}{2}x^2 + \frac{\alpha (\alpha -1)(\alpha -2)}{6}x^3 + ... + \left(\begin{array}{c}a\\ n\end{array}\right)x^n ", 
    strings_to_isolate = "="
)

f11 = MathTex(
    # ((alpha),(n)) = (alpha * (alpha - 1) * ... * (alpha -n+1))/ (n!)
    r"\left(\begin{array}{c}a\\ b\end{array}\right) = \frac{\alpha \cdot (\alpha -1) \cdot ... \cdot (\alpha -n +1)}{n!}", 
    strings_to_isolate = "="
)

f12 = MathTex(
    # ( 1 + x)^a = 1 + a*x
    r"(1+x)^\alpha \approx 1 + \alpha x", 
    substrings_to_isolate= "="
)


f13 = MathTex(
    # x = v^2/c^2 << 1
    r"x = \frac{\mathbf v^{2}}{c^{2}} \ll 1",
    substrings_to_isolate= "="
)


f14 = MathTex(
    # γ = (1 - v^2/c^2)^-1/2 
    r"\gamma = (1 - \frac{\mathbf v^2}{c^2})^{-\frac{1}{2}} ",
    # = 1 + -1/2* -v^2/c^2
    r"\approx 1 + (-\frac{1}{2})(- \frac{\mathbf v^2}{c^2})",
    # = 1 + v^2/c^2
    r" \approx 1 + \frac{\mathbf v^2}{c^2}",
    # & v<<c
    r" \wedge ( \mathbf v \ll c)",
    substrings_to_isolate= "="
)



f15 = MathTex(
    # Ek = mc^2(γ-1)
    r"E_{\kappa} = mc^2(\gamma -1)", 
    # = mc^2(1 + v^2/2c^2 -1)
    r" = mc^2(1+\frac{\mathbf v^2}{2c^2}-1)",
    # = 1/2 * m*v^2
    r"= \frac{1}{2}m\mathbf v^2",
    substrings_to_isolate= "="
)


class EnergiaCinetica(Scene):
    def construct(self):
        # f1=MathTex(r"{F} \cdot d\mathbf {x} "),                               # F*dx

        # f2 =MathTex(r"\mathbf {F} \cdot \mathbf {v} dt"),                     # F*v*dt

        # f3= MathTex(r"{\frac {d\mathbf {p} }{dt}}\cdot \mathbf {v} dt"),      # dp/dt*v*dt

        # f4 = MathTex(r"\mathbf {v} \cdot d\mathbf {p} "),                     # v*dp

        # f5 = MathTex(r"\mathbf {v} \cdot d(m\mathbf {v} )"),                  # v*d(m*v)

        # f6 = MathTex(r"= \text{d}( \mathbf {v} \cdot \mathbf {v} )"),         # d(v*v)

        # f7 = MathTex(r"=(d\mathbf {v} )\cdot \mathbf {v} +\mathbf {v} \cdot (d\mathbf {v} )"), # dv*v + v*dv

        # f8 = MathTex(r"=2(\mathbf {v} \cdot d\mathbf {v} )"),                 # 2(v*dv)

        # f9 = MathTex(r"={\displaystyle \mathbf {v} \cdot d(m\mathbf {v} )}"), # v*d(m*v)

        # f10 = MathTex(r"={\frac {m}{2}}d(\mathbf {v} \cdot \mathbf {v} )"),   # m/2*d(v*v)

        # f11 = MathTex(r"={\frac {m}{2}}dv^{2}"),                              # m/2*dv^2

        # f12 = MathTex(r"=d\left({\frac {mv^{2}}{2}}\right)"),                 # d(m*v^2/2)

        # F=m*a
        f13 = MathTex(r"{\displaystyle \mathbf {F}=m\mathbf {a} }"),

        # d(m*v)/dt
        f14 = MathTex(
            r" ={\frac {\mathrm {d} (m\mathbf {v} )}{\mathrm {d} t}}"),

        f15 = MathTex(
            r"=m\,{\frac {\,\mathrm {d} \mathbf {v} \,}{\mathrm {d} t}}"
        ),  # m*dv/dt

        f16 = MathTex(
            r"{\displaystyle \mathbf {a} = {\frac {\mathrm {d} (\mathbf {v} )}{\mathrm {d} t}}}"
        ),  # a = dv/dt

        # delta v / delta t
        f17 = MathTex(r"= \frac{\triangle \mathbf {v}}{\triangle t}"),

        # delta m = 0
        f18 = MathTex(r"\triangle m = 0"),

        f19 = MathTex(
            r"{\frac{\displaystyle \mathbf {F}}{\mathbf {a} }= m}")  # F / a =m

        # disegno F=m*a
        # sposto a sx
        fma = f13[0].center()
        self.play(Write(fma))
        self.play(fma.animate.to_edge(LEFT))

        # disegno  a dv/dt, delta v / delta t, m =0
        ta1 = MathTex(r"{\displaystyle \mathbf {a} = \mathrm {k}}").next_to(
            fma, RIGHT + UP)
        ta2 = f16[0].next_to(ta1, direction=RIGHT).shift(RIGHT)
        ta3 = f17[0].next_to(ta2, direction=RIGHT)
        tb1 = MathTex(r"{m} = \mathrm {k}").next_to(fma, RIGHT + DOWN)
        tb2 = f18[0].next_to(tb1,
                             direction=RIGHT).shift(np.array((0.8, 0.0, 0.0)))
        t = VGroup(ta1, tb1).shift(np.array((0.5, 0.0, 0.0)))

        # disegno graffa
        graffa = Brace(t, direction=LEFT, color=YELLOW, buff=0.05)

        self.play(GrowFromCenter(graffa))

        self.play(Write(ta1[0]))
        self.play(Write(ta2[0]))
        self.play(Write(ta3[0]))
        self.play(Write(tb1[0]))
        self.play(Write(tb2[0]))

        self.play(FadeOut(graffa), FadeOut(ta1[0]), FadeOut(ta2[0]),
                  FadeOut(ta3[0]), FadeOut(tb1[0]))
        self.play(fma.animate.to_edge(LEFT + UP), ta2.animate.to_edge(UP),
                  tb2.animate.to_edge(UP + RIGHT))

        # disegno 'L =' e
        # sposto a sinistra 'L'
        Ap_0_2 = f1[0:2].center()
        self.play(Write(Ap_0_2))
        self.play(Ap_0_2.animate.shift(LEFT))

        # disegno 'F*dx =' e lo sposto a sinistra
        # sposto a sx 'F ='
        Ap_2_4 = f1[2:4].next_to(Ap_0_2, direction=RIGHT)
        self.play(Write(Ap_2_4))
        self.play(Ap_0_2.animate.shift(LEFT), Ap_2_4.animate.shift(LEFT))

        # disegno 'F*v*dt' e evidenzio 'F'
        Ap_4_6 = f1[4:6].next_to(Ap_2_4, direction=RIGHT)
        self.play(Write(Ap_4_6))
        # disegno graffa sopra dx con scritto 's = v*t'
        graffa_spostamento = Brace(Ap_2_4[0][2:4],
                                   direction=UP,
                                   color=YELLOW,
                                   buff=0.05)
        definizone_spostamento = MathTex(
            r"\mathbf x=\mathbf {v} \cdot t").next_to(graffa_spostamento,
                                                      direction=UP)
        self.play(Indicate(Ap_2_4[0][0]), Indicate(Ap_4_6[0][0]))
        self.play(GrowFromCenter(graffa_spostamento),
                  GrowFromCenter(definizone_spostamento))
        # trasformo 'F = m* a' in 'F/a=m'
        # sposto a sxp_0_2,p_2_4,p_4_6, e elimino graffa_spostamento
        self.play(Ap_0_2.animate.shift(LEFT), Ap_2_4.animate.shift(LEFT),
                  Ap_4_6.animate.shift(LEFT),
                  FadeOutAndShift(graffa_spostamento, direction=LEFT),
                  FadeOutAndShift(definizone_spostamento, direction=LEFT))

        # disegno 'dp/dt*v*dt' e dopo crosso i due dt e elimino i deu dt e le cross e facccio replacetransform conp_8_10 ('v*dp')
        # sposto a sxp_0_2,p_2_4,p_4_6,p_8_10 (p_6_8 non c'e' perche ho fatto FadeOut)
        Ap_6_8 = f1[6:8].next_to(Ap_4_6, direction=RIGHT)
        self.play(Write(Ap_6_8))
        cross1 = Cross(Ap_6_8[0][3:5], stroke_width=4).scale(2)
        cross2 = Cross(Ap_6_8[0][8:10], stroke_width=4).scale(2)
        self.play(Create(cross1), Create(cross2))
        Ap_8_10 = f1[8:10].next_to(Ap_4_6, direction=RIGHT)
        self.play(ReplacementTransform(Ap_6_8.copy(), Ap_8_10),
                  FadeOut(Ap_6_8), ShrinkToCenter(cross1),
                  ShrinkToCenter(cross2))
        self.play(Ap_0_2.animate.shift(LEFT), Ap_2_4.animate.shift(LEFT),
                  Ap_4_6.animate.shift(LEFT), Ap_8_10.animate.shift(LEFT))
        # disegno graffa sopra dp con scritto 'p=m*v'
        graffa_quantita_moto = Brace(Ap_8_10[0:2],
                                     direction=UP,
                                     color=YELLOW,
                                     buff=0.05)
        definizone_quantita_moto = MathTex(
            r"\mathbf {p}= m \cdot \mathbf {v} ").next_to(graffa_quantita_moto,
                                                          direction=UP)
        self.play(GrowFromCenter(graffa_quantita_moto),
                  GrowFromCenter(definizone_quantita_moto))

        # disegno 'v*d(m*v)' NB: f1[10:11] è 11 e non 12 perche ho tolto un meno alla fine della formual (che nelle altre invece c'è)
        # riduco la fomula a 'L = v*d(m*v)'
        Ap_10_12 = f1[10:11].next_to(Ap_8_10, direction=RIGHT)
        self.play(Write(Ap_10_12))
        # creao VGroup con 'L = v*d(p)'
        gruppo_formula_per_relativita = VGroup(Ap_0_2[0][0:1], Ap_8_10)
        formula_per_relativita = f5[0:3].to_edge(DOWN + LEFT).set_color(PINK)
        self.play(
            ReplacementTransform(gruppo_formula_per_relativita.copy(),
                                 formula_per_relativita))
        # sposto a sx Ap_0_2, Ap_10_12 e elimino graffa_qunatita_di_moto
        self.remove(Ap_8_10)
        self.play(
            FadeOut(Ap_2_4), FadeOutToPoint(Ap_4_6, LEFT),
            Ap_10_12.animate(run_time=1.5).next_to(Ap_0_2, direction=RIGHT),
            FadeOutAndShift(graffa_quantita_moto, direction=LEFT),
            FadeOutAndShift(definizone_quantita_moto, direction=LEFT))

        # disegno 'v(dm*v + m*dv)'
        # sposto a sx Ap_0_2, Ap_10_12, Bp_0_2 e elimino graffa_delta_m
        Bp_0_2 = f2[0:2].next_to(Ap_10_12, direction=RIGHT)
        self.play(ReplacementTransform(Ap_10_12.copy(), Bp_0_2), run_time=1.5)
        # disegno graffa sopra dm con scritto 'delta m / delta t= 0'
        graffa_delta_m = Brace(Bp_0_2[1][2:4],
                               direction=UP,
                               color=YELLOW,
                               buff=0.05)
        definizone_delta_m = MathTex(
            r"\frac{\triangle m}{\triangle t}=0").next_to(graffa_delta_m,
                                                          direction=UP)
        self.play(GrowFromCenter(graffa_delta_m),
                  GrowFromCenter(definizone_delta_m))
        self.play(Indicate(tb2[0], color=RED_E),
                  Indicate(definizone_delta_m[0][0:2], color=RED_E))
        self.play(Ap_0_2.animate.shift(LEFT), Ap_10_12.animate.shift(LEFT),
                  Bp_0_2.animate.shift(LEFT),
                  FadeOutAndShift(graffa_delta_m, direction=LEFT),
                  FadeOutAndShift(definizone_delta_m, direction=LEFT))

        # disegno 'v*m*dv' NB: f2[2:3] è 3 e non 4 perche ho tolto un meno alla fine della formual (che nelle altre invece c'è)
        # sposto a sx Ap_0_2, Ap_10_12, Bp_0_2, Bp_2_4
        Bp_2_4 = f2[2:4].next_to(Bp_0_2, direction=RIGHT)
        self.play(ReplacementTransform(Bp_0_2.copy(), Bp_2_4), run_time=1.5)
        self.play(FadeOut(Ap_10_12),
                  FadeOutToPoint(Bp_0_2, LEFT),
                  Bp_2_4.animate.next_to(Ap_0_2[0:1], direction=RIGHT),
                  FadeOut(Ap_0_2[1:2]),
                  run_time=1.5)

        # sposto dL = m*v*dv in alto a sinistra
        L1 = VGroup(Ap_0_2[0:1], Bp_2_4)
        self.add(L1)
        self.clear()
        self.add(formula_per_relativita.to_edge(DOWN + LEFT))
        self.play(L1.animate.move_to(3.3 * LEFT + 2.8 * UP))

        # disegno 'd(v*v)' e sposto a sx
        Cp_0_2 = f3[0:2].center()
        self.play(Write(Cp_0_2))
        self.play(Cp_0_2.animate.shift(LEFT))

        # disegno 'dv*v + v*dv'
        # sposto a sx Cp_0_2 e 'dv*v + v*dv' (Cp_2_4)
        Cp_2_4 = f3[2:4].next_to(Cp_0_2, direction=RIGHT)
        self.play(ReplacementTransform(Cp_0_2.copy(), Cp_2_4))
        self.play(Cp_0_2.animate.shift(LEFT), Cp_2_4.animate.shift(LEFT))

        # disegno '2(v*dv)'
        # sposto a sx Cp_0_2, Cp_2_4, Cp_4_6
        Cp_4_6 = f3[4:6].next_to(Cp_2_4, direction=RIGHT)
        self.play(ReplacementTransform(Cp_2_4.copy(), Cp_4_6))
        self.play(Cp_0_2.animate.shift(LEFT), Cp_2_4.animate.shift(LEFT),
                  Cp_4_6.animate.shift(LEFT))

        # elimino 'dv*v + v*dv' e riduco formula a 'd(v*v) = 2(v*dv)'
        self.play(
            FadeOutToPoint(Cp_2_4, LEFT),
            Cp_4_6.animate(run_time=1.5).next_to(Cp_0_2, direction=RIGHT))
        L2 = VGroup(Cp_0_2, Cp_4_6)
        self.add(L2)
        self.remove(Cp_0_2, Cp_4_6)
        self.play(L2.animate.move_to(3.3 * RIGHT + 2.8 * UP))

        # divido per 1/2 d(v*v)=2(v*dv)
        L3 = MathTex(
            r"\frac{1}{2}\text{d}( \mathbf {v} \cdot \mathbf {v} ) =(\mathbf {v} \cdot d\mathbf {v} ) "
        ).next_to(L2, aligned_edge=LEFT).shift(2 * LEFT)
        self.play(ReplacementTransform(L2, L3))

        # # disegno rettangoli attorno a v*dv
        rettangolo1 = SurroundingRectangle(Bp_2_4[1][2:6],
                                           color=YELLOW,
                                           buff=0.20)
        rettangolo2 = SurroundingRectangle(Cp_4_6[10:16],
                                           color=YELLOW,
                                           buff=0.20)
        self.play(Create(rettangolo1), run_time=1.5)
        self.play(Create(rettangolo2), run_time=1.5)

        # disegno L = m/2d(v^2) = Ek
        Dp_0_2 = f4[0][0:2].center()
        self.play(Write(Dp_0_2))
        self.play(Dp_0_2.animate.shift(LEFT))
        Dp_2_4 = f4[1:3].next_to(Dp_0_2, direction=RIGHT)
        self.play(Write(Dp_2_4))
        self.play(Dp_0_2.animate.shift(LEFT), Dp_2_4.animate.shift(LEFT))

        # traslo'dL=d(m/2*v^2)' in basso a dx
        formula_lavoro = VGroup(Dp_0_2, Dp_2_4)
        self.play(formula_lavoro.animate.to_edge(DOWN + RIGHT).set_color(TEAL))

        # diseno Ek =
        # sposto Ek a sx
        Ep_0_2 = f4a[0:2].center()
        self.play(Write(Ep_0_2))
        self.play(Ep_0_2.animate.shift(LEFT * 1.4))

        # disegno ∫F*dx
        # sposto a sx Ek, ∫F*dx
        Ep_2_4 = f4a[2:4].next_to(Ep_0_2, direction=RIGHT)
        self.play(Write(Ep_2_4))
        self.play(Ep_0_2.animate.shift(LEFT * 1.4),
                  Ep_2_4.animate.shift(LEFT * 1.4))

        # disegno '∫v*d(m*v)'
        # sposto a sx Ek, ∫F*dx, ∫v*d(m*v)
        Ep_4_6 = f4a[4:6].next_to(Ep_2_4, direction=RIGHT)
        self.play(Write(Ep_4_6))
        self.play(Ep_0_2.animate.shift(LEFT * 1.4),
                  Ep_2_4.animate.shift(LEFT * 1.4),
                  Ep_4_6.animate.shift(LEFT * 1.4))

        # disegno '∫d(m/2*v^2)'
        # sposto a sx Ek, ∫F*dx, ∫v*d(m*v), ∫d(m/2*v^2)
        Ep_6_8 = f4a[6:8].next_to(Ep_4_6, direction=RIGHT)
        self.play(Write(Ep_6_8))
        self.play(Ep_0_2.animate.shift(LEFT * 1.4),
                  Ep_2_4.animate.shift(LEFT * 1.4),
                  Ep_4_6.animate.shift(LEFT * 1.4),
                  Ep_6_8.animate.shift(LEFT * 1.4))

        # disegno 'm*v^2/2'
        # sposto a sx Ek, ∫F*dx, ∫v*d(m*v), ∫d(m/2*v^2), m*v^2/2
        Ep_8_10 = f4a[8:10].next_to(Ep_6_8, direction=RIGHT)
        self.play(Write(Ep_8_10))
        # evidenzio '∫d(m/2*v^2)' con formula_lavoro
        self.play(Indicate(Ep_6_8, color=TEAL),
                  Indicate(formula_lavoro, color=TEAL))
        # evidenzio 'Ek' e 'm/2*v^2'
        self.play(Indicate(Ep_0_2[0:1], color=RED_E),
                  Indicate(Ep_8_10, color=RED_E))
        # riduco fomrula f4a(gruppo_formula_energia_cinetica_classica) a 'Ek = m/2*v^2 = ∫L'(formula_energia_cinetica)
        gruppo_formula_energia_cinetica_classica = VGroup(
            Ep_0_2, Ep_2_4, Ep_4_6, Ep_6_8, Ep_8_10)
        formula_energia_cinetica_classica = MathTex(
            r"E_{\kappa} = {\frac {mv^{2}}{2}} =\int d\!L").center()
        self.play(
            ReplacementTransform(gruppo_formula_energia_cinetica_classica,
                                 formula_energia_cinetica_classica))
        # sposto i alto a sx formula_energia_cinetica
        # elimino gruppo_formula_energia_cinetica_classica, elimino L1, L2, elimino i 2 rettangoli, elimino formula_lavoro
        self.play(FadeOutAndShift(rettangolo1, direction=UP),
                  FadeOutAndShift(rettangolo2, direction=UP),
                  FadeOutAndShift(L1, direction=UP),
                  FadeOutAndShift(L3, direction=UP),
                  FadeOutAndShift(formula_lavoro, direction=DOWN),
                  formula_energia_cinetica_classica.animate.to_edge(UP + LEFT))

        # disegno 'E = ∫L ='
        # sposto a sx
        Fp_0_2 = f6[0:2].center()  # Ek =
        Fp_2_4 = f6[2:4].next_to(Fp_0_2, direction=RIGHT)  # ∫L =
        self.play(Write(Fp_0_2))
        self.play(Write(Fp_2_4))
        self.play(Fp_0_2.animate.shift(LEFT), Fp_2_4.animate.shift(LEFT))
        # evidenzio '∫L' in formula_energia_cinetica_classica e in Fp_2_4
        self.play(Indicate(formula_per_relativita, color=YELLOW),
                  Indicate(Fp_2_4, color=YELLOW))
        # disegno '∫(v*dp) =' e elimino Fp_2_4 e copio formula_per_relativita in Fp_4_6
        Fp_4_6 = f6[4:6].next_to(Fp_0_2)
        self.play(FadeOut(Fp_2_4),
                  ReplacementTransform(formula_per_relativita.copy(), Fp_4_6))

        # disegno in alto a destra 'p=m*γ*v'
        quantita_di_moto = p.to_edge(UP + RIGHT)
        self.play(Write(quantita_di_moto))

        # disegno fattore lorenziano in basso a dx
        fattore_lorenziano = gamma.to_edge(DOWN + RIGHT)
        self.play(Write(fattore_lorenziano))

        # disegno '∫v*(m*γ*v)'
        Fp_6_8 = f6[6:8].next_to(Fp_4_6, direction=RIGHT)
        # evidenzio quantita_di_moto e Fp_4_6
        self.play(Indicate(quantita_di_moto, color=ORANGE),
                  Indicate(Fp_4_6[0][4], color=ORANGE))
        # trasporto quantita_di_moto in Fp_6_8
        self.play(ReplacementTransform(quantita_di_moto.copy(), Fp_6_8))
        # spost a sx Fp_0_2, Fp_4_6, Fp_6_8 NB: NON C'E' Fp_2_4 perche ho fatto il Transformreplace la riga sopra
        self.play(Fp_0_2.animate.shift(5 * LEFT),
                  Fp_4_6.animate.shift(5 * LEFT),
                  Fp_6_8.animate.shift(5 * LEFT))
        # disegno 'm*γ*v - ∫m*γ*v * dv'
        Fp_8_10 = f6[8:10].next_to(Fp_6_8)
        self.play(Write(Fp_8_10))
        # trasformo Fp_8_10 in Fp_10_12
        Fp_10_12 = f6[10:12].next_to(Fp_6_8)
        self.play(ReplacementTransform(Fp_8_10, Fp_10_12))

        # elimino '∫(v*dp) =', '∫v*(m*γ*v)'
        self.play(FadeOutAndShift(Fp_4_6, LEFT), FadeOutAndShift(Fp_6_8, LEFT),
                  Fp_10_12.animate.next_to(Fp_0_2))

        # evidenzio il fattore_lorenziano e γ in Fp_10_12
        self.play(Indicate(fattore_lorenziano, color=YELLOW),
                  Indicate(Fp_10_12[0][1], color=YELLOW, scale_factor=2),
                  Indicate(Fp_10_12[0][9], color=YELLOW, scale_factor=2))

        # trasformo Fp_10_12 in Fp_12_14('m*γ*v^2 - (-m*c^2/2)*∫γ*d(1 - v^2/c^2) =')
        Fp_12_14 = f6[12:14].next_to(Fp_10_12)
        self.play(ReplacementTransform(Fp_10_12.copy(), Fp_12_14))

        # elimino Fp_10_12
        # sposto a sx Fp_12_14
        self.play(FadeOutAndShift(Fp_10_12, LEFT),
                  Fp_12_14.animate.next_to(Fp_0_2))

        # trasformo 'm*γ*v^2 - (-m*c^2/2)*∫γ*d(1 - v^2/c^2) ='(Fp_12_14) in 'm*γ*v^2 + m*c^2*(1 - v^2/c^2)^1/2 - E0' (Fp_14_16)
        Fp_14_16 = f6[14:16].next_to(Fp_0_2)
        self.play(ReplacementTransform(Fp_12_14, Fp_14_16))

        # trasformo 'm*γ*v^2 + m*c^2*(1 - v^2/c^2)^1/2 - E0' (Fp_14_16) in 'm*γ*( v^2 + c^2*(1 - v^2/c^2)) - E0' (Fp_16_18)
        Fp_16_18 = f6[16:18].next_to(Fp_0_2)
        self.play(ReplacementTransform(Fp_14_16, Fp_16_18))

        # trasformo 'm*γ*( v^2 + c^2*(1 - v^2/c^2)) - E0' (Fp_16_18) in 'm*γ*( v^2 + c^2 - v^2 ) - E0' (Fp_18_20)
        Fp_18_20 = f6[18:20].next_to(Fp_0_2)
        self.play(ReplacementTransform(Fp_16_18, Fp_18_20))

        # trasformo 'm*γ*( v^2 + c^2 - v^2 ) - E0' (Fp_18_20) in 'm*γ*c^2 - E0' (Fp_20_22)
        Fp_20_22 = f6[20:22].next_to(Fp_0_2)
        self.play(ReplacementTransform(Fp_18_20, Fp_20_22))

        # disegno '(γ -1)* m*c^2' Fp_22_24
        Fp_22_24 = f6[22:24].next_to(Fp_20_22)
        self.play(Write(Fp_22_24))

        # evidenzio Fp_0_2('E') e Fp_22_24 ('(γ -1)* m*c^2')
        self.play(Indicate(Fp_0_2, color=RED_E), Indicate(Fp_22_24,
                                                          color=RED_E))

        # creo vgroup con fomrula finale  e formula finale con cui fare trasformazione
        formula_energia_cinetica_relativistica = MathTex(
            r"E_{\kappa} = (\gamma -1) mc^{2} ").center()
        gruppo_formula_energia_cinetica_relativistica = VGroup(
            Fp_0_2, Fp_20_22, Fp_22_24)

        # trasformo gruppo_formula_energia_cinetica_relativistica in formula_energia_cinetica_relativistica
        self.play(
            ReplacementTransform(gruppo_formula_energia_cinetica_relativistica,
                                 formula_energia_cinetica_relativistica))

        self.wait(2)
        # disegno 'm*γ*v^2 - (-m*c^2/2)*∫γ*d(1 - v^2/c^2)'




class GraficoInversaClassica(GraphScene):
    def __init__(self, **kwargs):
        GraphScene.__init__(self,
                            y_min=0,
                            x_min=0,
                            y_max= 3.5,
                            x_max= 3.5,
                            graph_origin= 3 * DOWN + 4 * LEFT,#    5 * DOWN + 6 * LEFT,
                            y_axis_label= r"$E_{\kappa} \rightleftharpoons {\mathbf v} $",
                            x_axis_label= r"${\mathbf v}  \rightleftharpoons E_{\kappa}$",
                            y_axis_config={
                                "tick_frequency": 7,
                                "include_tip": True,
                                # "numbers_with_elongated_ticks": [1],
                                # "numbers_to_exclude": [1,2,3]
                            },
                            x_axis_config={
                                "tick_frequency": 7,
                                "include_tip": True,
                                # "numbers_to_exclude": [1,2,3]
                                
                            },
                            **kwargs)
    
    
    def construct(self):
       
        self.setup_axes(animate=True)
        m = 1
        velocita_luce = 5
        
        
        # DISEGNO GRAFICO ENERGIA CLASSICA E SUA FUNZIONE INVERSA + LINNE ORIZZONTALI PER VERIFICARE CHE E' BIUNIVOCA
        grafico_energia_classica = self.get_graph(lambda x : 0.5*m*math.pow(x,2), color = BLUE, x_min=0, x_max = 2.7)
        grafico_energia_classica.set_stroke(width=5)
        label_grafico_energia_classica = self.get_graph_label(grafico_energia_classica,
                                       x_val=3,
                                       color=BLUE,
                                       direction=LEFT,
                                       buff=SMALL_BUFF,
                                       label=r"{E_{\kappa}} = \frac{m \mathbf v^{2}}{2}").shift(LEFT*1.5)
        
        grafico_energia_classica_invertita = self.get_graph(lambda x : math.sqrt(2*x/m), x_min =0, x_max = 3.3, color = ORANGE)
        grafico_energia_classica_invertita.set_stroke(width=5)
        label_grafico_energia_classica_invertita = self.get_graph_label(grafico_energia_classica_invertita,
                                       x_val=3,
                                       color=ORANGE,
                                       direction=DOWN ,
                                       buff=MED_LARGE_BUFF,
                                       label=r"\mathbf {v} = \sqrt{\frac{2 \cdot {E_{\kappa}}}{m}}")
        
        grafico_bisettrice_primo_terzo_quadrante = self.get_graph(lambda x :x, x_min =0, color = GREY, x_max = 3.1)
        grafico_bisettrice_primo_terzo_quadrante.set_stroke(width=3)
        
        
        # disegno grafico Ek e label
        self.play(Create(grafico_energia_classica, lag_ratio=0.03, run_time=1.5), Create(label_grafico_energia_classica, lag_ratio=0.03, run_time=1.5))
        
        # disegno ed elimino linee per verificare che funzinoe sia biunivoca
        giu = DOWN + UP *4 
        linee = VGroup()
                
        for i in range(5):
            giu += DOWN
            linee.add(DashedLine(config.left_side + RIGHT*3.1+ giu , config.right_side +  LEFT*2 + giu  ).set_color(RED_E))
        
        for j in range(5):
            self.play(Create(linee[j]))
            
        for k in range(5):
            self.play(Uncreate(linee[k]))
            
        
        #disegno Ek invertita + label e  bisettrice
        self.play(Create(grafico_energia_classica_invertita, lag_ratio=0.03, run_time=1.5), Create(label_grafico_energia_classica_invertita, lag_ratio=0.03, run_time=1.5))
        self.play(Create(grafico_bisettrice_primo_terzo_quadrante, lag_ratio=0.03, run_time=1.5))
        
        # trasformo la invertitaa in quella normale e viceversa 
        self.play(ClockwiseTransform(grafico_energia_classica_invertita.copy(), grafico_energia_classica))
        self.wait(0.5)
        self.play(CounterclockwiseTransform(grafico_energia_classica.copy(), grafico_energia_classica_invertita))
        

        
        
        
        
        
        


class GraficoInversaRelativistica(GraphScene):
    
    
    
    def __init__(self, **kwargs):
        GraphScene.__init__(self,
                            y_min=0,
                            x_min=0,
                            y_max=20,
                            x_max=21,
                            graph_origin= 3 * DOWN + 4 * LEFT,
                            y_axis_label= r"$E_{\kappa} \rightleftharpoons {\mathbf v} $",
                            x_axis_label= r"${\mathbf v}  \rightleftharpoons E_{\kappa}$",
                            y_axis_config={
                                "tick_frequency": 40,
                                "include_tip": True
                            },
                            x_axis_config={
                                "tick_frequency": 42,
                                "numbers_with_elongated_ticks": [c],
                                "include_tip": True,
                                
                            },
                            **kwargs)

    def construct(self):
    
        # DISEGNO GRAFICO ENERGIA RELATIVISTICA E SUA FUNZIONE INVERSA + LINNE ORIZZONTALI PER VERIFICARE CHE E' BIUNIVOCA
        
        self.setup_axes(animate=True)
        
        c = 20
        m = 0.005
        
        grafico_energia_relativistica = self.get_graph(lambda x:  ((1 / (math.sqrt(1 - math.pow( x, 2) / math.pow(c, 2)))) - 1) * m * math.pow(c, 2), color = RED,  x_min=0, x_max=c - 0.1,)
        
        grafico_energia_relativistica.set_stroke(width=5)
        label_grafico_energia_relativistica = self.get_graph_label(grafico_energia_relativistica,
                                       x_val=10,
                                       color=RED,
                                       direction= UP+RIGHT,
                                       buff=SMALL_BUFF,
                                       label=r"{E_{\kappa}} = \frac{mc^{2}}{\sqrt{1-\frac{\mathbf v^{2}}{c^{2}}}}-mc^{2}").shift(LEFT*1.5).scale(0.7)
        
        grafico_energia_relativistica_invertita = self.get_graph(lambda x:   c*math.sqrt(math.pow(x,2)/ (math.pow(x,2) + math.pow(m,2)*math.pow(c,4))), x_min =0, x_max = c-1, color = GREEN)
        grafico_energia_relativistica_invertita.set_stroke(width=5)
        label_grafico_energia_relativistica_invertita = self.get_graph_label(grafico_energia_relativistica_invertita,
                                       x_val=6,
                                       color=GREEN,
                                       direction=DOWN ,
                                       buff=MED_LARGE_BUFF,
                                       label=r"{\mathbf v} = c\sqrt{\frac{E_{\kappa}^2}{E_{\kappa}^2+m^2c^4}}  ").scale(0.7)
        
        grafico_bisettrice_primo_terzo_quadrante = self.get_graph(lambda x :0.97*x, color = GREY, x_min = 0, x_max = c-1)
        grafico_bisettrice_primo_terzo_quadrante.set_stroke(width=3)
        
        
        # disegno grafico Ek e label
        self.play(Create(grafico_energia_relativistica, lag_ratio=0.03, run_time=1.5), Create(label_grafico_energia_relativistica, lag_ratio=0.03, run_time=1.5))
        
        
        # disegno ed elimino linee per verificare che funzinoe sia biunivoca
        giu = DOWN + UP *4 
        linee = VGroup()
                
        for i in range(5):
            giu += DOWN
            linee.add(DashedLine(config.left_side + RIGHT*3.1+ giu , config.right_side +  LEFT*2 + giu  ).set_color(YELLOW))
        
        for j in range(5):
            if(j==4):
                self.play(FadeOut(label_grafico_energia_relativistica, run_time= 0.8), Create(linee[j]))
                continue
            self.play(Create(linee[j]))
            
        for k in range(5):
            if(k==4):
                self.play(Create(label_grafico_energia_relativistica, run_time= 1.5), Uncreate(linee[k]))
                continue
            self.play(Uncreate(linee[k]))
            
            
            
        #disegno Ek invertita + label e  bisettrice e asindoto in 'c'
        self.play(Create(grafico_energia_relativistica_invertita, lag_ratio=0.03, run_time=1.5), Create(label_grafico_energia_relativistica_invertita, lag_ratio=0.03, run_time=1.5))
        
        c_line = self.get_vertical_line_to_graph(c - 0.1,
                                                 grafico_energia_relativistica,
                                                 line_class=DashedLine,
                                                 color=RED)
        c_letter = Text('c', color=RED).next_to(c_line, direction=DOWN)
        
        self.play(Create(c_line), Write(c_letter))
        
        self.play(Create(grafico_bisettrice_primo_terzo_quadrante))
        
        
        # trasformo la invertitaa in quella normale e viceversa 
        self.play(ClockwiseTransform(grafico_energia_relativistica_invertita.copy(), grafico_energia_relativistica))
        self.wait(0.5)
        self.play(CounterclockwiseTransform(grafico_energia_relativistica.copy(), grafico_energia_relativistica_invertita))
        


class FunzioniInverseConLimiti( GraphScene, MovingCameraScene):
    
    def __init__(self, **kwargs):
        GraphScene.__init__(self,
                            y_min=0,
                            x_min=0,
                            y_max= 3.5,
                            x_max= 7,
                            graph_origin= 5 * DOWN + 6 * LEFT,
                            y_axis_label= r"${\mathbf v} $",
                            x_axis_label= r"$E_{\kappa}$",
                            y_axis_config={
                                "tick_frequency": 17,
                                "include_tip": True,
                                "numbers_with_elongated_ticks": [1],
                                "numbers_to_exclude": [1,2,3]
                            },
                            x_axis_config={
                                "tick_frequency": 7,
                                "include_tip": True,
                                "numbers_to_exclude": [1,2,3]
                                
                            },
                            **kwargs)
    
    
    def construct(self):
        

        # FUNZIONI INVERSA ENERGIA CLASSICA E LIMITE
        G_p_0_1 = f7[0:1].center()  # Ek = m/2*v^2
        G_p_1_2 = f7[1:2].center()  # 2Ek=m/2*v^2
        G_p_2_3 = f7[2:3].center()  # 2Ek/m = v^2
        G_p_3_4 = f7[3:4].center()  # v = √(2*ek)/m
        G_p_4_5 = f7[4:5].center()  # lim Ek->∞ di v(Ek) =  √(2EK/m)
        # lim Ek->∞ di v(Ek) =  √(2(+ ∞)/m)
        G_p_5_6 = f7[5:6].center().shift(0.3 * RIGHT)
        G_p_6_7 = f7[6:7].center().shift(0.3 * DOWN + 0.3 * LEFT)  # lim Ek->∞ di v(Ek) -> + ∞
        self.play(ReplacementTransform(G_p_0_1, G_p_1_2))
        self.remove(G_p_0_1)
        self.play(ReplacementTransform(G_p_1_2, G_p_2_3))
        self.remove(G_p_1_2)
        self.play(ReplacementTransform(G_p_2_3, G_p_3_4))
        self.remove(G_p_2_3)
        # sposto in alto a sx v = √(2*ek)/m (G_p_3_4)
        self.play(G_p_3_4.animate.to_edge(UP + LEFT))
        # trasformo (copiando) v = √(2*ek)/m (G_p_3_4) in lim Ek->∞ di v(Ek) =  √(2EK/m) (G_p_4_5)
        self.play(ReplacementTransform(G_p_3_4.copy(), G_p_4_5))
        # trasforo formule a aprtire da dopo '=' del limite
        self.play(ReplacementTransform(G_p_4_5[0][11:], G_p_5_6[0][11:]))
        # trasforo formule a aprtire da dopo '=' del limite
        self.play(ReplacementTransform(G_p_5_6[0][11:], G_p_6_7[0][11:]))
        
        # elimino ttuo dallo schermo e trasloin alto a dx lim Ek->∞ di v(Ek) -> + ∞
        self.clear()
        self.add(G_p_3_4.to_edge(UP + LEFT), G_p_6_7)
        self.play(G_p_6_7.animate.to_edge(UP + RIGHT))
        
        
        # FUNZIONI INVERSA ENERGIA RELATIVISTICA E LIMITE
        H_p_0_1 = f8[0:1].center() # Ek = (γ -1)* m*c^2
        H_p_1_2 = f8[1:2].center() # mc^2/√(1-v^2/c^2) -mc^2
        H_p_2_3 = f8[2:3].center() # Ek + mc^2 = mc^2*γ
        H_p_3_4 = f8[3:4].center() # Ek/mc^2 + 1 = γ
        H_p_4_5 = f8[4:5].center() # Ek/mc^2 + 1 = 1/√(1-v^2/c^2)
        H_p_5_6 = f8[5:6].center() # [(Ek+mc^2)/mc^2]^-1 = √(1-v^2/c^2)
        H_p_6_7 = f8[6:7].center() # mc^2/(Ek = mc^2) = √(1-v^2/c^2)
        H_p_7_8 = f8[7:8].center() # - m^2*c^4/(Ek^2 + m^2*c^4) +1 = v^2/c^2
        H_p_8_9 = f8[8:9].center() # [-m^2c^4 + Ek^2 + m^2c^4]/[Ek^2 +m^2c^4] = v^2/c^2
        H_p_9_10 = f8[9:10].center() # Ek^2/[Ek^2 +m^2c^4] = v^2/c^2
        H_p_10_11 = f8[10:11].center() # c^2 * Ek^2/[Ek^2 +m^2c^4] = v^2
        H_p_11_12 = f8[11:12].center() # c √(Ek^2/[Ek^2 +m^2c^4]) = v
        H_p_12_13 = f8[12:13].center() # lim Ek->∞ di v(Ek) =  c^2 √(Ek^2/[Ek^2 +m^2c^4])
        H_p_13_14 = f8[13:14].center().shift(RIGHT*0.2+DOWN*0.1) # c √ {Ek^2/Ek^2}* {1+[1+m^2c^4/Ek^2]}
        H_p_14_15 = f8[14:15].center().shift(LEFT*0.2+DOWN*0.1) # c √ {Ek^2/Ek^2}* {1+[1+m^2c^4/∞}
        H_p_15_16 = f8[15:16].center().shift(LEFT*0.6) # c √ {Ek^2/Ek^2}* {1+[1+0}
        H_p_16_17 = f8[16:17].center().shift(LEFT*0.2+DOWN*0.1) # c =  3*10^8 m/s
        
        #scrivo 'mc^2/√(1-v^2/c^2) -mc^2' e la trasformo fino a farla diventare 'c √(Ek^2/[Ek^2 +m^2c^4]) = v'
        self.play(Write(H_p_0_1))
        self.play(ReplacementTransform(H_p_0_1, H_p_1_2))
        self.play(ReplacementTransform(H_p_1_2, H_p_2_3))
        self.play(ReplacementTransform(H_p_2_3, H_p_3_4))
        self.play(ReplacementTransform(H_p_3_4, H_p_4_5))
        self.play(ReplacementTransform(H_p_4_5, H_p_5_6))
        self.play(ReplacementTransform(H_p_5_6, H_p_6_7))
        self.play(ReplacementTransform(H_p_6_7, H_p_7_8))
        self.play(ReplacementTransform(H_p_7_8, H_p_8_9))
        self.play(ReplacementTransform(H_p_8_9, H_p_9_10))
        self.play(ReplacementTransform(H_p_9_10, H_p_10_11))
        self.play(ReplacementTransform(H_p_10_11, H_p_11_12))
        
        # sposto 'c √(Ek^2/[Ek^2 +m^2c^4]) = v' in basso a sx
        self.play(H_p_11_12.animate.to_edge(LEFT + DOWN))
    
        # scrivo al centro 'lim Ek->∞ di v(Ek) =  c^2 √(Ek^2/[Ek^2 +m^2c^4])'
        self.play(Write(H_p_12_13[0][0:12].copy()))        
        self.play(Write(H_p_12_13[0][12:]))
        
        # trasformo 'lim Ek->∞ di v(Ek) =  c^2 √(Ek^2/[Ek^2 +m^2c^4])' in 'c =  3*10^8 m/s' (NB: trsaformo solo la parte dopo "lim Ek-> ∞ di v(Ek) =")
        self.play(ReplacementTransform(H_p_12_13[0][12:], H_p_13_14[0][12:]))
        self.play(ReplacementTransform(H_p_13_14[0][12:], H_p_14_15[0][12:]))
        self.play(ReplacementTransform(H_p_14_15[0][12:], H_p_15_16[0][12:]))
        self.play(ReplacementTransform(H_p_15_16[0][12:], H_p_16_17[0][12:]))

        # elimino tutto dallo schermo
        self.clear()
        
        
        # ridisegno tutto come prima del 'self.clear()' : basso-sx 'c √(Ek^2/[Ek^2 +m^2c^4]) = v' - basso-dx 'lim Ek->∞ di v(Ek) = c =  3*10^8 m/s' - alto-sx 'v = √(2*ek)/m' - alto-dx 'lim Ek->∞ di v(Ek) -> + ∞'
        self.add(G_p_3_4.to_edge(UP + LEFT), G_p_6_7.to_edge(UP + RIGHT),H_p_11_12.to_edge(LEFT + DOWN), H_p_16_17)

        # sposto in basso a dx 'lim Ek->∞ di v(Ek) = c =  3*10^8 m/s'
        self.play(H_p_16_17.animate.to_edge(DOWN+RIGHT))
        
        # elimino 'v = √(2*ek)/m' & 'c √(Ek^2/[Ek^2 +m^2c^4]) = v' + sposto in alto-sx lim Ek->∞ di v(Ek) -> + ∞ (G_p_6_7) e in alto-dx c =  3*10^8 m/s (H_p_16_17)
        self.play(FadeOutAndShift(G_p_3_4, direction = LEFT), FadeOutAndShift(H_p_11_12, direction = LEFT), G_p_6_7.animate.to_edge(UP+LEFT), H_p_16_17.animate.to_edge(UP+RIGHT))
        
        # rimpicciolisco il frame e sposto in alto i limiti 
        self.play(self.camera.frame.animate.set_width(20), G_p_6_7.animate.shift(UP).set_color(BLUE), H_p_16_17.animate.shift(UP).set_color(YELLOW) )
        
        
        self.setup_axes(animate=True)
        m = 1
        velocita_luce = 1
        
        grafico_energia_classica = self.get_graph(lambda x : math.sqrt(2*x/m), color = BLUE, x_min=0)
        grafico_energia_classica.set_stroke(width=5)
        label_grafico_energia_classica = self.get_graph_label(grafico_energia_classica,
                                       x_val=3,
                                       color=BLUE,
                                       direction=UP + LEFT,
                                       buff=SMALL_BUFF,
                                       label=r"\mathbf {v} = \sqrt{\frac{2 \cdot {E_{\kappa}}}{m}}")
        
        grafico_energia_relativistica = self.get_graph(lambda x : velocita_luce**2*math.sqrt(math.pow(x,2)/ (math.pow(x,2) + math.pow(m,2)*math.pow(velocita_luce,4))), x_min =0, color = YELLOW)
        grafico_energia_relativistica.set_stroke(width=5)
        label_grafico_energia_relativistica = self.get_graph_label(grafico_energia_relativistica,
                                       x_val=5,
                                       color=YELLOW,
                                       direction=UP ,
                                       buff=MED_LARGE_BUFF,
                                       label=r"{\mathbf v} = c\sqrt{\frac{E_{\kappa}^2}{E_{\kappa}^2+m^2c^4}}  ")
        
        self.play(Create(grafico_energia_classica, lag_ratio=0.03, run_time=1.5), Create(label_grafico_energia_classica, lag_ratio=0.03, run_time=1.5))
        self.play(Create(grafico_energia_relativistica, lag_ratio=0.03, run_time=1.5), Create(label_grafico_energia_relativistica, lag_ratio=0.03, run_time=1.5))
        
        giu = DOWN*3.3
        c_line = DashedLine(config.left_side + RIGHT*1.4 + giu, config.right_side +  LEFT*4.1 + giu ).set_color(RED_E)
        label_c_line = Text("c").next_to(c_line, direction = LEFT).shift(LEFT*0.5).set_color(RED_E)
        self.play(Create(label_c_line), Create(c_line))
        

        
       
        
        

class ConfrontoGrafici(GraphScene, MovingCameraScene):
    
    def __init__(self, **kwargs):
        GraphScene.__init__(self,
                            y_min=0,
                            x_min=0,
                            y_max=c / 9,
                            x_max=c + 10,
                            y_axis_label="$m$",
                            x_axis_label="$v$",
                            y_axis_config={
                                "tick_frequency": 1000, #c / 100,
                                "include_tip": True
                            },
                            x_axis_config={
                                "tick_frequency": 1000, #c / 10,
                                # "numbers_with_elongated_ticks": [c],
                                "include_tip": True,
                                
                            },
                            **kwargs)

    def construct(self):
        m = 0.005
        
        self.setup_axes(animate=True)
        self.camera.frame.save_state()

        # create reference graphs for the moving dots
        retta_riferimento = self.get_graph(lambda x: 0.07 * x,
                                           color=BLACK,
                                           x_min=0,
                                           x_max=c + 10)

        KRel_riferimento = self.get_graph(lambda x: 0.6 * ((1 / (math.sqrt(
            1 - math.pow(x, 2) / math.pow(c, 2)))) - 1) * m * math.pow(c, 2),
                                          x_min=0,
                                          x_max=c - 1,
                                          color=BLACK)
        # create displayed graphs

        K = self.get_graph(lambda x: 0.3 * m * math.pow(x, 2),
                           color=BLUE,
                           x_min=0,
                           x_max=c + 10)

        KRel = self.get_graph(lambda x: 0.6 * ((1 / (math.sqrt(1 - math.pow(
            x, 2) / math.pow(c, 2)))) - 1) * m * math.pow(c, 2),
                              x_min=0,
                              x_max=c - 0.1,
                              color=RED_E)

        # add labels to graphs
        K_label = self.get_graph_label(K,
                                       x_val=c + 5,
                                       color=BLUE,
                                       direction=DOWN + RIGHT,
                                       buff=SMALL_BUFF,
                                       label=r"{E_{\kappa}}=\frac{1}{2}mv^{2}")
        KRel_label = self.get_graph_label(
            KRel,
            x_val=c / 2 + 3,
            color=RED_E,
            buff=SMALL_BUFF,
            direction=LEFT + UP,
            label=r"{E_{\kappa}}=mc^{2}(\gamma-1)")

        # create dots
        moving_dot_k = Dot().move_to(K.points[0]).set_color(YELLOW)
        moving_dot_KRel = Dot().move_to(
            KRel_riferimento.points[0]).set_color(YELLOW)
        moving_dot_retta_riferimento = Dot(radius=0.000000001).move_to(
            retta_riferimento.points[0]).set_color(BLACK)

        self.add(retta_riferimento)
        self.play(Write(K, lag_ratio=0.03, run_time=1.5),
                  Write(K_label, lag_ratio=0.03, run_time=1.5))
        self.play(Write(KRel, lag_ratio=0.03, run_time=1.5),
                  Write(KRel_label, lag_ratio=0.03, run_time=1.5))

        c_line = self.get_vertical_line_to_graph(c - 0.1,
                                                 KRel,
                                                 line_class=DashedLine,
                                                 colore=RED_E)
        c_letter = Text('c', color=RED_E).next_to(c_line, direction=DOWN)
        self.play(Create(c_line), Write(c_letter))

        self.play(self.camera.frame.animate.scale(0.7).move_to(moving_dot_k))

        def update_curve(mob):
            mob.move_to(moving_dot_retta_riferimento.get_center())

        self.camera.frame.add_updater(update_curve)
        self.play(
            MoveAlongPath(moving_dot_retta_riferimento,
                          retta_riferimento,
                          rate_func=linear,
                          run_time=4),
            MoveAlongPath(moving_dot_k, K, rate_func=linear, run_time=4.5),
            MoveAlongPath(moving_dot_KRel,
                          KRel_riferimento,
                          rate_func=linear,
                          run_time=7))
        self.camera.frame.remove_updater(update_curve)

        self.play(Restore(self.camera.frame))
        
       
       
        
       


class RiemannRectangles(GraphScene):
    def __init__(self, **kwargs):
        GraphScene.__init__(self,
                            # y_min=0,
                            # x_min=0,
                            y_max= 8,
                            x_max= 6,
                            y_axis_height= 6,
                            graph_origin= 2.5 * DOWN + 5 * LEFT,#    5 * DOWN + 6 * LEFT,
                            y_axis_label= r"$ F$",
                            x_axis_label= r"$ s$",
                            y_axis_config={
                                # "tick_frequency": 7,
                                # "include_tip": True,
                                # "numbers_with_elongated_ticks": [1],
                                # "numbers_to_exclude": [1,2,3]
                                "include_numbers": False,
                            },
                            x_axis_config={
                                "tick_frequency": 0.5,
                                # "include_tip": True,
                                # "numbers_to_exclude": [1,2,3]4
                                "include_numbers": False,
                                "numbers_with_elongated_ticks": [1, 5],
                                
                            },
                            **kwargs)
    
    
    def construct(self):
       
        self.setup_axes(animate=True)
        
        graph = self.get_graph(lambda x : 0.05 *(x+3-5) *(x+3-5)* (x-5) *(x+5)+5  , color = YELLOW, x_min=0.2, x_max = 5.5) #  0.3 * (x+ 3-5) * (x+ 3-5) * (x-5) + 5
        graph.set_stroke(width=5)
        label_graph = self.get_graph_label(graph,
                                       x_val=3,
                                       color=YELLOW,
                                       direction=UP+LEFT,
                                       buff=SMALL_BUFF,
                                       label=r"d \!L = f(s)").shift(UP*0.8 +LEFT)
        
        self.play(Create(graph), Create(label_graph), run_time = 1.5)
        
        def rect(x):
            return x
        
        recta = self.get_graph(rect, x_min = -1, x_max = 5)
        
        kwargs = {
            "x_min" : 1,
            "x_max" : 5, 
            "fill_opacity" : 0.75, 
            "stroke_width" : 0.25,
            
        }
        
        self.graph = graph 
        iterazioni = 6
        
        self.rect_list = self.get_riemann_rectangles_list(graph, iterazioni, start_color= BLUE, end_color=ORANGE, **kwargs)
            
        
        
        
        a = Text("a").next_to(self.get_vertical_line_to_graph(1, graph), direction = DOWN).scale(0.8).set_color(BLUE)
        b = Text("b").next_to(self.get_vertical_line_to_graph(5, graph), direction = DOWN).scale(0.8).set_color(ORANGE)
        
        
        integrale = MathTex("\int_{a}^{b} f(s)ds").center().shift(RIGHT * 4 + UP ).scale(1.2)
        
        integrale[0][2].set_color(BLUE)
        integrale[0][1].set_color(ORANGE)
        
        self.play(Write(integrale[0][0]))
        self.play(Write(integrale[0][2]))
        self.play(Write(integrale[0][1]))
        
        self.play(ReplacementTransform(integrale[0][2].copy(), a ))
        self.play(ReplacementTransform(integrale[0][1].copy(), b ))
        
        
        for rect in self.rect_list[0]:
                    self.play(GrowFromEdge(rect, edge = DOWN))
        
        tacca = VGroup(self.get_vertical_line_to_graph(2, graph), self.get_vertical_line_to_graph(2.5, graph))
        
        graffa_dx = Brace(tacca, direction = DOWN)    
        
        dx = MathTex("ds").next_to(graffa_dx, direction = DOWN)
        
   
        
        altezza = self.get_vertical_line_to_graph(2, graph)
        
        
        fs_brace = Brace(altezza, LEFT, buff =0 )
        
        fs_label = fs_brace.get_text("$f(s)$", buff = SMALL_BUFF)
        
        
        self.play(
            self.rect_list[0].animate.set_fill(opacity =0.3),
            self.rect_list[0][2].animate.set_fill(opacity = 1),
            GrowFromCenter(fs_brace),
            Write(fs_label)
        )
        
        self.play(Write(graffa_dx), Write(dx))
        
        self.play(ReplacementTransform(fs_label.copy(), integrale[0][3:7]))
        self.play(ReplacementTransform(dx.copy(), integrale[0][7:9].shift(RIGHT*0.1)))
            
        
        scale_factor = [0.25, 0.35, 0.37, 0.39, 0.41, 0.43]
        
        for j in range(1,6):                
            
            self.transform_between_riemann_rects(self.rect_list[j-1], self.rect_list[j], replace_mobject_with_target_in_scene = True, dx = 1, run_time= 0.9)
            self.add(fs_brace, fs_label)
            
        
            graffa_dx_loop = VGroup()
            
            graffa_dx_loop.add(Brace(VGroup(self.get_vertical_line_to_graph(2, graph), self.get_vertical_line_to_graph(2.5 - scale_factor[j-1], graph)), direction = DOWN))
        
            dx_loop = VGroup()
            
            dx_loop.add(MathTex("ds").next_to(dx, direction = OUT).scale(1 - j *0.1))
            
            
            if (j==1):
                self.play(Transform(graffa_dx, graffa_dx_loop[0]), Transform(dx, MathTex("dx").next_to(dx, direction = OUT).scale(1 - j *0.1)))    
            else:
                self.play(Transform( graffa_dx, Brace(VGroup(self.get_vertical_line_to_graph(2, graph), self.get_vertical_line_to_graph(2.5 - scale_factor[j-1], graph)), direction = DOWN)), Transform(dx, MathTex("dx").next_to(graffa_dx, direction = OUT).scale(1 - j *0.1).shift(DOWN*0.5)))
                
                
        ds_tendente = Brace(integrale, direction= DOWN)
        ds_label = ds_tendente.get_text("$ ds \\rightarrow 0$").next_to(ds_tendente, direction = DOWN)
        
        self.play(GrowFromCenter(ds_tendente), Write(ds_label))
           
        L = MathTex(r" =\!L").next_to(integrale).shift(LEFT*0.5)
        area = MathTex(r"\!L").center().shift(LEFT*1.1 + DOWN * 1.1 ).scale(1.5)
        
        self.play(FadeOut(ds_tendente), FadeOut(ds_label), integrale.animate.shift(LEFT*0.5) ,Write(L))
        
        self.play( FadeOut(fs_brace), FadeOut(fs_label), FadeIn(area))
          
        self.play(ApplyWave(VGroup(integrale, L, area), color = RED), ApplyWave(self.rect_list[-1]), run_time = 2)
        
            
        
    






class TaylorMcLaurin(GraphScene, MovingCameraScene):
    def __init__(self, **kwargs):
        GraphScene.__init__(self,
                            y_min=-1,
                            x_min=-3,
                            y_max= 3,
                            x_max=3,
                            graph_origin = ORIGIN + DOWN * 5 + LEFT * 1.1,
                            # y_axis_label= r"${\mathbf v} $",
                            x_axis_label= r"$x$",
                            y_axis_config={
                                "tick_frequency": 1,
                                # "include_tip": True,
                            },
                            x_axis_config={
                                "tick_frequency": 1,
                                # "include_tip": True,
                                "include_numbers": True,
                                "numbers_to_show": [1],
                                "numbers_with_elongated_ticks": [1],
                                # "numbers_to_exclude": [2,3]
                                
                            },
                            **kwargs)
    
    
    def construct(self):
        print(self.camera.frame.get_width())
        self.camera.frame.save_state()
        self.camera.frame.set_width(25)

        self.setup_axes(animate=True)
        
        
        grafico_originale = self.get_graph(lambda x : math.sqrt(1+x), color = BLUE, x_min = -1, x_max = 3)
        label_grafico_originale = self.get_graph_label(grafico_originale, x_val = 2, color = BLUE, direction = UP + LEFT, buff = SMALL_BUFF, label = r"\sqrt{1+x}")
        
        self.play(Write(grafico_originale, lag_ratio=0.03, run_time=1.5), Write(label_grafico_originale, lag_ratio=0.03, run_time=1.5))
        
        
        # disegno f9 pezzo per pezzo 
        taylor_cinetica = f9.center().shift(UP * 2)
        self.play(ReplacementTransform(label_grafico_originale, taylor_cinetica[0:1].set_color(BLUE))) # √1+x
        self.play(Write(taylor_cinetica[1:2])) # (1 + x )^1/2 
        self.play(Write(taylor_cinetica[2:3]))
        self.play(Write(taylor_cinetica[3:4]))
        self.play(Write(taylor_cinetica[4:5].set_color(YELLOW))) # 1
        self.play(Write(taylor_cinetica[5:6].set_color(RED)))    # + x/2
        self.play(Write(taylor_cinetica[6:7].set_color(PURPLE))) # - x^2/8
        self.play(Write(taylor_cinetica[7:8].set_color(ORANGE))) # + x^3/16
        self.play(Write(taylor_cinetica[8:9].set_color(PINK)))   # - 5*x^4/128
        
        taylor_generale = f10.next_to(taylor_cinetica.get_center(), direction = UP).shift(UP+ RIGHT * 4)
        coefficiente_binomiale_generalizzato = f11.next_to(taylor_generale, direction = UP).set_color(GREEN_SCREEN)
        
        
        self.play(ReplacementTransform(taylor_cinetica[1:9].copy(), taylor_generale))
        self.play(Indicate(taylor_generale[0][7], color= YELLOW), Indicate(taylor_cinetica[4:5], color= YELLOW))
        self.play(Indicate(taylor_generale[0][8:11], color= RED), Indicate(taylor_cinetica[5:6], color= RED))
        self.play(Indicate(taylor_generale[0][11:22], color= PURPLE), Indicate(taylor_cinetica[6:7], color= PURPLE))
        self.play(Indicate(taylor_generale[0][22:38], color= ORANGE), Indicate(taylor_cinetica[7:8], color= ORANGE))
        
        self.play(Indicate(taylor_generale[0][43:], color = GREEN_SCREEN), ApplyMethod(taylor_generale[0][43:].set_color, GREEN_SCREEN))
        self.play(ReplacementTransform(taylor_generale[0][43:].copy(), coefficiente_binomiale_generalizzato))
        self.play(Write(Text("Sviluppo di Taylor-Mc Laurin del binomio").next_to(taylor_generale, direction = LEFT).scale(0.7).shift(RIGHT*1.8)))
        
        self.play(self.camera.frame.animate.shift(2*DOWN) )
        self.play(self.camera.frame.animate.set_width(17))
        
        grafici = VGroup()
        graffe = VGroup()
        colori = [YELLOW, RED, PURPLE, ORANGE, PINK]
        
        funzioni = [( lambda x: 1), (lambda x: 1+ x/2), (lambda x: 1 + x/2- math.pow(x, 2)/8),(lambda x: 1 + x/2- math.pow(x, 2)/8 + math.pow(x, 3)/16),(lambda x: 1 + x/2- math.pow(x, 2)/8 + math.pow(x, 3)/16 - 5*math.pow(x, 4)/128)]
        
        for i in range(5):
            if (i==0):
                  graffe.add(Brace(taylor_cinetica[4:5+i],direction=DOWN,   color=WHITE )).scale(1.5)
            else:
                graffe.add(Brace(taylor_cinetica[4:5+i],direction=DOWN,   color=WHITE ))
            
            grafici.add(self.get_graph(funzioni[i], color = colori[i]))
            
            if (i==0):
                self.play(Create(grafici[i]), GrowFromCenter(graffe[i]))
            else:
                self.play(ReplacementTransform(graffe[i-1], graffe[i]))
                self.play(ReplacementTransform(grafici[i-1], grafici[i]))
                
                
        self.play(self.camera.frame.animate.shift(15*RIGHT), taylor_cinetica.animate.shift(RIGHT*15))         
        
        sviluppo_potenza = f12.next_to(taylor_cinetica[2], direction = DOWN).shift(RIGHT*1.1 + DOWN*0.8)
        
        self.play(ReplacementTransform(taylor_cinetica[2].copy(), sviluppo_potenza))
        
        x_come_rapport_fra_velocita = f13.next_to(f12, direction = RIGHT).shift(UP*0.05 + RIGHT*0.3)
        
        self.play(ReplacementTransform(sviluppo_potenza[0][-1].copy(), x_come_rapport_fra_velocita))
        
        self.play(CreateThenDestruction(MathTex(r"\frac{\mathbf v^{2}}{c^{2}} = \frac{37,5}{3 \cdot 299792458} = 0,000000125").next_to(taylor_cinetica, direction = DOWN).shift(DOWN*4)))
        
        gamma_binomio = f14.next_to(taylor_cinetica, direction = DOWN).shift(DOWN*2)
        
        self.play(ReplacementTransform(sviluppo_potenza[0][0:6].copy(), gamma_binomio[0:3]))
        self.play(ReplacementTransform(gamma_binomio[2:3].copy(), gamma_binomio[3:4]))
        self.play(ReplacementTransform(gamma_binomio[3:4].copy(), gamma_binomio[4:5]))
        
        self.play(Write(gamma_binomio[5:6]))
        
        formula_con_approssimazinoe = f15.next_to(gamma_binomio, direction = DOWN).shift(DOWN)
        
        self.play(Write(formula_con_approssimazinoe[0:3]))
        
        self.play(Indicate(gamma_binomio[0]), Indicate(gamma_binomio[4]), Indicate(formula_con_approssimazinoe[0:3][2][4]), run_time = 1.5)
        
        self.play(ReplacementTransform(formula_con_approssimazinoe[0:3].copy(), formula_con_approssimazinoe[3:6]))
        
        cross1 = Cross(formula_con_approssimazinoe[3:6][2][4]).scale(1.8)
        cross2 = Cross(formula_con_approssimazinoe[3:6][2][13]).scale(1.8)
        self.play(Create(cross1), Create(cross2))
        self.play(ShrinkToCenter(cross1), ShrinkToCenter(cross2), ReplacementTransform(formula_con_approssimazinoe[3:6], formula_con_approssimazinoe[6:9].next_to(formula_con_approssimazinoe[0:3])))
               
        
        
        


def functioncurlreal(p, velocity=1):
    x, y = p[:2]
    result =  - y * RIGHT + x * UP
    result *= velocity
    return result


class VectorField(Scene):
    def construct(self):
        # func = lambda pos: np.sin(pos[1]/2)*RIGHT+np.cos(pos[0]/2)*UP
        vector_field = ArrowVectorField(functioncurlreal)
        vector_field.set_color(BLUE)
        # self.add(vector_field)
        self.play(*[GrowArrow(vec) for vec in vector_field])
        # circle = Circle(radius=2).shift(LEFT)
        # self.add(circle.copy().set_color(GRAY))
        dot = Dot([1,1.5,0], color=ORANGE).scale(3)
        dot.add(Text("+").move_to(dot.get_center()))
        
        circle = Circle(radius=1.9)
        self.play(Create(circle))
        
        # vector_field.nudge(circle, -2, 60, True)
        vector_field.nudge(dot, -2, 60)

        # circle.add_updater(vector_field.get_nudge_updater(pointwise=True))
        dot.add_updater(vector_field.get_nudge_updater())
        self.add( dot)
        self.wait(1)
        
        



class Titolo(Scene):
    def construct(self):
        titolo = Text("Energia Cinetica")
        self.play(Write(titolo), run_time= 1.5)
        self.play(FadeOut(titolo), run_time = 1.5)

        
class Titolo1(Scene):
    def construct(self):
        titolo = Text("Confronto Grafici")
        self.play(Write(titolo), run_time= 1.5)
        self.play(FadeOut(titolo), run_time = 1.5)

class Titolo2(Scene):
    def construct(self):
        titolo = Text("Funzioni Inverse & Limiti")
        self.play(Write(titolo), run_time= 1.5)
        self.play(FadeOut(titolo), run_time = 1.5)
        
class Titolo3(Scene):
    def construct(self):
        titolo = Text("Grafico Funzione Inversa (Classica)")
        self.play(Write(titolo), run_time= 1.5)
        self.play(FadeOut(titolo), run_time = 1.5)

class Titolo4(Scene):
    def construct(self):
        titolo =Text("Grafico Funzione Inversa (Relativistica)")
        self.play(Write(titolo), run_time= 1.5)
        self.play(FadeOut(titolo), run_time = 1.5)
        
        
class Titolo5(Scene):
    def construct(self):
        titolo = Text("Cos'è un Integrale?")
        self.play(Write(titolo), run_time= 1.5)
        self.play(FadeOut(titolo), run_time = 1.5)
        

class Titolo6(Scene):
    def construct(self):
        titolo = Text("Taylor - Mc Laurin")
        self.play(Write(titolo), run_time= 1.5)
        self.play(FadeOut(titolo), run_time = 1.5)
        
        
        
class Titolo7(Scene):
    def construct(self):
        titolo = Text("Particella in campo magnetico")
        self.play(Write(titolo), run_time= 1.5)
        self.play(FadeOut(titolo), run_time = 1.5)
        
        
        
class Titolo8(Scene):
    def construct(self):
        titolo = Text("Effetti relativistici")
        self.play(Write(titolo), run_time= 1.5)
        self.play(FadeOut(titolo), run_time = 1.5)
        
        


class MotoParticella(Scene):
    def construct(self):
        formula1 = MathTex(r"F = ma = qE")
        formula2 = MathTex(r"\mathbf v = \frac{2 \pi r}{T}")
        formula3 = MathTex(r"r = \frac{m \mathbf v}{qB} = \frac{m}{qB}\cdot \frac{2 \pi r}{T}")
        formula4 = MathTex(r"T =\frac{2 \pi m}{qB}")
        self.play(Write(formula1.center().shift(UP*2+LEFT*3)))
        self.play(Write(formula2.center().shift(UP*2+RIGHT*3)))
        self.play(Write(formula3.center()))
        self.play(Write(formula4.center().shift(DOWN*2)))
        
        
class MotoParticellaRelativistica(Scene):
    def construct(self):
        formula1 = MathTex(r"F = \frac{dp}{dt}= \frac{d(\gamma (t) \cdot m \cdot \mathbf v(t))}{dt}= m \mathbf v \cdot \frac{d \gamma}{dt} + m \gamma \cdot\frac{d \mathbf v}{dt}")
        formula2 = MathTex(r"\frac{d \gamma}{dt}=[\frac{1}{\sqrt{1-\frac{\mathbf v^{2}}{c^2}}}]' = - \frac{\gamma \mathbf v a}{c^2 - \mathbf v^2}")
        formula3 = MathTex(r"F = m \mathbf v  a(1 - \frac{ \mathbf v^2 c^2}{c^2 -  \mathbf v^2} )= m \gamma a \cdot \frac{c^2 -  \mathbf v^2 c^2}{c^2 -  \mathbf v^2} = q  \mathbf v B")
        formula4 = MathTex(r"a = \frac{(q \mathbf v B) \cdot (c^2 - \mathbf v^2)}{m \gamma c^2 ( 1 - \mathbf v^2)}")
        self.play(Write(formula1.center().shift(UP *3)))
        self.play(Write(formula2.center().shift(UP*0.8)))
        self.play(Write(formula3.center().shift(DOWN)))
        self.play(Write(formula4.center().shift(DOWN * 3)))
        
        
        






#  for(i=0; i<x;i++){ // 2x+2
#      for(j=o; j<x+i; j++){ // x(x+1)/2 
#         esegui istruzioni... 
#      }
#  }