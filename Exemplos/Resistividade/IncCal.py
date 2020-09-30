from sympy import *
import numpy as np
init_printing()


def get_vars(expression):
  return expression.free_symbols

def val_real(expression,pares_val,prec):
  return round(expression.evalf(subs=pares_val),prec)

def get_dict_expr_inct(expression):
  vars = get_vars(expression)
  pares_val = dict()
  pares_inc = dict()
  unids = dict()
  for var in vars:
    aux = float(input("Digite o valor da variável "+ str(var)+": "))
    pares_val[var] = aux
    aux1 = float(input("E sua respectiva incerteza: "))
    pares_inc[var] = aux1
    aux2 = str(input("Agora as suas unidades: "))
    unids[var] = aux2
  return pares_val, pares_inc,unids
 
def incert_relativa(expression,pares_val,pares_inc,prec):
  vars = get_vars(expression)
  valor = 0
  for var in vars:
    deriv = float(expression.diff(var).evalf(subs = pares_val))
    val = pares_inc[var]
    value = (val * deriv)**2
    valor += value
  resp = round(np.sqrt(valor),prec)
  return resp

def der_latex(expression,vari,pares_val,pares_inc,i,unidades,unids_vars,file):
    vars = get_vars(expression)
    file.write("Os valores das variáveis são:\n \\begin{align*} \n")
    count = 0
    for var in vars:
        if count == 0:
            file.write(str(latex(sympify(var)))+" &= ("+str(pares_val[var])+"\\pm "+str(pares_inc[var])+") "+ unids_vars[var])
        if count != 0:
            file.write("\\\\ "+str(latex(sympify(var)))+" &= ("+str(pares_val[var])+"\\pm "+str(pares_inc[var])+") "+ unids_vars[var])
        count += 1
    file.write("\n \\end{align*} \n")
        
    file.write("As derivadas parcias da função "+str(vari)+"("+str(vars)+"), substituindo os valores das variáveis, são: \n")
    file.write("\\begin{align}"+"\n")
    valores = dict()
    valf = round(expression.evalf(subs=pares_val),i)
    count = 0
    for var in vars:
        deriv = str(latex(diff(expression,var)))
        if(count != 0):
            stringExp = "\\\\ \\frac{\partial "+vari+"}{\\partial "+str(latex(sympify(var)))+"} = " + deriv
        if(count == 0):
            stringExp = " \\frac{\partial "+vari+"}{\\partial "+str(latex(sympify(var)))+"} = " + deriv
        file.write(stringExp+ "&= ")
        a = round(diff(expression,var).evalf(subs = pares_val),i)
        valores[var] = a
        der = str(latex(a))
        file.write(der)
        count += 1
    file.write("\n"+"\\end{align}"+"\n")
    file.write("Calculando agora a incerteza, tem-se: \n")
    file.write("\\begin{align} \n \\Delta "+str(vari)+" &= \\sqrt{")
    count = 0
    for var in vars:
        if(count != 0):
            file.write(" + \\left( \\frac{\\partial "+str(vari)+"}{\\partial "+str(latex(sympify(var)))+"}\\Delta "+str(latex(sympify(var)))+"\\right)^{ 2 }")
        if(count == 0):
            file.write("\\left( \\frac{\\partial "+str(vari)+"}{\\partial "+str(latex(sympify(var)))+"}\\Delta "+str(latex(sympify(var)))+"\\right)^{ 2 }")
        count+=1
    file.write("} \\\\ \\Delta "+str(vari)+" &= \\sqrt{")
    count = 0
    re = 0
    for var in vars:
        re += (valores[var] * pares_inc[var]) ** 2
        if(count != 0):
            file.write(" + \\left("+str(valores[var])+"\\times "+ str(pares_inc[var])+"\\right)^{ 2 }")
        if(count == 0):
            file.write("\\left("+str(valores[var])+"\\times "+ str(pares_inc[var])+"\\right)^{ 2 }")
        count+=1
    file.write("} \\\\ \\Delta "+str(vari)+" &= \\sqrt{")
    count = 0
    file.write(str(latex(re)))
    r = sqrt(re)
    r = round(r,i)
    file.write("} \\\\ \\Delta "+str(vari)+" &= "+str(r)+" "+str(units) + "\n \\end{align}")
    valf,r,mult = ajuste(valf,r,i)
    file.write("\n Portanto, o valor final de "+str(vari)+ " é:\n \\begin{equation} \n ")
    if mult == 0:
        file.write("\\boxed{ "+str(vari)+" = \left("+str(valf)+" \\pm "+str(r)+"\\right) "+str(units)+"}\n \\end{equation}")
    if mult != 0 and isinstance(mult,str) == false:
        file.write("\\boxed{ "+str(vari)+" = \left("+str(valf)+" \\pm "+str(r)+"\\right) 10^{ "+str(mult)+"} "+str(units)+"}\n \\end{equation}")
    if mult != 0 and isinstance(mult,str):
        file.write("\\boxed{ "+str(vari)+" = \left("+str(valf)+" \\pm "+str(r)+"\\right) "+str(mult)+str(units)+"}\n \\end{equation}")

prefixos = {-6:"\\mu ",-3:"m ",3:"k ", 6:"M "}



def ajuste(resul,incert,ida):
    for i in range(-6,6):
        if((incert / (10**i)) <= 10 and (incert / (10**i)) >= 1):
            if (i!=1) and (i != 0) and (i!=-1):
                guess = i
                nincert = round((incert/(10**i)),1)
                nresul = round(resul / (10 ** i),1)
            else: 
                guess = 0
                nresul = resul
                nincert = incert
    if guess in (-6,-3,3,6):
        guess = prefixos[guess]
    return nresul,nincert,guess
    
    


def inicio_docu(expression,vari,file):
    section = " \\section{Cálculo da incerteza de " + str(vari)+ "}"
    intro = "Para o cálculo da incerteza de "+str(vari)+", usaremos o método das \\emph{Derivadas Parcias}. A seguinte equação descreve a função "+str(vari)+": \\\\"
    beg = "\\begin{equation} \n "+str(vari)+" = " + str(latex(expression)) + "\n" + "\\end{equation} \n" 
    file.write(section)
    file.write(intro)
    file.write(beg)

nome = str(input("Digite o nome do Experimento: "))
f = open(nome+".tex","w")

vari = str(input("Digite aqui a variável a ser calculada: "))
expr = sympify(str(input("Digite aqui sua expressão para {0}: ".format(vari))))
units = str(input("Digite aqui a unidade da variável {0}: ".format(vari)))
expr_lat = str(latex(expr))
parvals, parinc,unidads = get_dict_expr_inct(expr)

i = int(input("Digite a precisão dos cálculos: "))

f.write("\\documentclass{article}" + "\n"+"\\usepackage{amsmath}"+"\n"+"\\begin{document}")
inicio_docu(expr,vari,f)
der_latex(expr,vari,parvals,parinc,i,units,unidads,f)
f.write("\\end{document}")
f.close()
