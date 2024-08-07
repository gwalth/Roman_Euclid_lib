
def original(xpix, ypix, xl, yl):
    box_hw_x = int(xl/2)
    box_hw_y = int(yl/2)

    if xpix > 0:
        xpix = int(xpix+0.5)
    else:
        xpix = int(xpix-0.5)

    if ypix > 0:
        ypix = int(ypix+0.5)
    else: 
        ypix = int(ypix-0.5)

    y0 = ypix-box_hw_y
    y1 = ypix+box_hw_y
    x0 = xpix-box_hw_x
    x1 = xpix+box_hw_x

    print(x0, x1, y0, y1, xl, yl, x1-x0, y1-y0)


def new(xpix, ypix, xl, yl):

    if xpix > 0:
        xsign = 1
    else:
        xsign = 1  

    if ypix > 0:
        ysign = 1
    else:
        ysign = 1  

    y0 = int(ypix - yl/2. + ysign * 0.5)
    y1 = int(ypix + yl/2. + ysign * 0.5)
    x0 = int(xpix - xl/2. + xsign * 0.5)
    x1 = int(xpix + xl/2. + xsign * 0.5)

    print(x0, x1, y0, y1, xl, yl, x1-x0, y1-y0)

def new_test(xpix, ypix, xl, yl):
    if xpix > 0:
        xsign = 1
    else:
        xsign = -1  

    if ypix > 0:
        ysign = 1
    else:
        ysign = -1 


    #print(xpix, ypix, xpix_new, ypix_new)

    # either 0.0 or 0.5 fractions of pixel
    if xl % 2:
        # if remainder
        xpix_new = int(xpix) + xsign*0.5

    else:
        # if no remainder
        xpix_new = int(xpix + xsign*0.5)

    if yl % 2:
        # if remainder
        ypix_new = int(ypix) + ysign*0.5

    else:
        # if no remainder
        ypix_new = int(ypix + ysign*0.5)

    box_hw_x = xl/2
    box_hw_y = yl/2

    y0 = int(ypix_new - box_hw_y)
    y1 = int(ypix_new + box_hw_y)
    x0 = int(xpix_new - box_hw_x)
    x1 = int(xpix_new + box_hw_x)

    dx = xpix - xpix_new
    dy = ypix - ypix_new

    print(x0, x1, y0, y1, xl, yl, x1-x0, y1-y0, dx, dy)
    #print(xpix, ypix, xpix_new, ypix_new)




    #y0 = int(ypix - yl/2. + ysign * 0.5)
    #y1 = int(ypix + yl/2. + ysign * 0.5)
    #x0 = int(xpix - xl/2. + xsign * 0.5)
    #x1 = int(xpix + xl/2. + xsign * 0.5)

    #print(x0, x1, y0, y1, xl, yl, x1-x0, y1-y0)


def kevin_test(xpix, ypix, xl, yl):

    x0 = int(round(xpix-0.5)+1.0-xl/2.)
    x1 = x0 + xl
    y0 = int(round(ypix-0.5)+1.0-yl/2.)
    y1 = y0 + yl

    xpix_new = (x0 + x1)/2.
    ypix_new = (y0 + y1)/2.

    dx = xpix - xpix_new
    dy = ypix - ypix_new

    print(x0, x1, y0, y1, xl, yl, x1-x0, y1-y0, dx, dy)

####################################

#original(1000.0, 1000.000, 21, 21)

#new(1000.0, 1000.000, 21, 21)

#new_test(1000.0, 1000.000, 21, 21)

#kevin_test(1000.0, 1000.000, 21, 21)


xpix = [1000.0, 1000.0,   1000.0,   -1.0, -1.0,   -1.0,   -1.0,   -1.0,   -1.0,   -1.0  ]
ypix = [1000.0, 1000.499, 1000.999, -1.0, -1.499, -1.999, -1.000, -1.222, -1.499, -1.999]
xl   = [21,     21,       21,        20,   20,     20,     21,     21,     21,     21]
yl   = [21,     21,       21,        20,   20,     20,     21,     21,     21,     21]


for i in range(len(xpix)):

    print(     xpix[i], ypix[i], xl[i], yl[i])
    original(  xpix[i], ypix[i], xl[i], yl[i]) # fails
    new(       xpix[i], ypix[i], xl[i], yl[i]) # fails
    new_test(  xpix[i], ypix[i], xl[i], yl[i])
    kevin_test(xpix[i], ypix[i], xl[i], yl[i])
    print()





