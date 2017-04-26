# the basic reaction equations
def meinhardt2_1_a(oldA, oldI, p):
    dA = p['sdens_0'] * (((oldA * oldA) / oldI) + p['proda']) - p['rema'] * oldA
    return dA

def meinhardt2_1_i(oldA, oldI, p):
    dI = p['sdens_0'] * oldA * oldA - p['remi'] * oldI + p['prodi']
    return dI

#an activator-substrate systems
def meinhardt_paper1_a(oldA, oldI, p):
    aStar2 = ((oldA*oldA)/(1 + p['sata'] * oldA*oldA)) + p['proda']
    dA = (p['sdens'] * oldI * aStar2) - p['rema'] * oldA
    return dA

def meinhardt_paper1_i(oldA, oldI, p):
    aStar2 = ((oldA*oldA)/(1 + p['sata'] * oldA*oldA)) + p['proda']
    dI = p['prodi'] - p['sdens'] * oldI * aStar2 - p['remi'] * oldI
    return dI

#an activator-inhibitor system
def meinhardt_paper2_a(oldA, oldI, p):
    aStar2 = (oldA*oldA)/(1 + p['sata'] * oldA*oldA)
    dA = (p['sdens'] * (aStar2 + p['proda']))/oldI - p['rema'] * oldA
    return dA

def meinhardt_paper2_i(oldA, oldI, p):
    aStar2 = (oldA*oldA)/(1 + p['sata'] * oldA*oldA);
    dI = p['sdens'] * aStar2 - p['remi'] * oldI + p['prodi'];
    return dI

# a dual inhibitor system
def meinhardt5_2_a(oldA, oldI, oldI2, p):
    aStar2 = (oldA*oldA + p['proda']) / (1 + p['sata'])
    dA = p['sdens']*aStar2*(oldI+oldI2) - p['rema']*oldA
    return dA

def meinhardt5_2_i(oldA, oldI, oldI2, p):
    aStar2 = (oldA*oldA + p['proda']) / (1 + p['sata'])
    dI = p['prodi'] - p['sdens']*oldI*aStar2 - p['remi']*oldI
    return dI

def meinhardt5_2_i2(oldA, oldI, oldI2, p):
    aStar2 = (oldA*oldA + p['proda']) / (1 + p['sata'])
    dI2 = p['prodi2'] - p['sdens']*oldI2*aStar2 - p['remi2']*oldI2
    return dI2

# a dual inhibitor system
def meinhardt5_3_a(oldA, oldI, oldI2, p):
    dA = (p['sdens']/oldI2) * (oldA*oldA/oldI + p['proda']) - p['rema']*oldA
    return dA

def meinhardt5_3_i(oldA, oldI, oldI2, p):
    dI = p['remi']*oldA*oldA/oldI2 - p['remi']*oldI + p['prodi']
    return dI

def meinhardt5_3_i2(oldA, oldI, oldI2, p):
    dI2 = p['remi2']*(oldA - oldI2)
    return dI2

# an alternate dual inhibitor system
def meinhardt5_4_a(oldA, oldI, oldI2, p):
    dA = p['sdens'] * (oldA*oldA + p['proda'])/(p['sati']*oldI + p['sati2']*oldI2) - p['rema'] * oldA
    return dA

def meinhardt5_4_i(oldA, oldI, oldI2, p):
    dI = p['remi'] * (oldA*oldA - oldI) + p['prodi']
    return dI

def meinhardt5_4_i2(oldA, oldI, oldI2, p):
    dI2 = p['remi2'] * (oldA - oldI2)
    return dI2
