
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
####################################

###original(1000.0, 1000.000, 21, 21)
#original(1000.0, 1000.499, 21, 21)
#original(1000.0, 1000.999, 21, 21)
#original(-1.0, -1.000, 20, 20)
#original(-1.0, -1.499, 20, 20)
#original(-1.0, -1.999, 20, 20)
#original(-1.0, -1.000, 21, 21)
#original(-1.0, -1.499, 21, 21)
#original(-1.0, -1.999, 21, 21)

#print()

#new(1000.0, 1000.000, 21, 21)
#new(1000.0, 1000.499, 21, 21)
#new(1000.0, 1000.999, 21, 21)
#new(-1.0, -1.000, 20, 20)
#new(-1.0, -1.499, 20, 20)
#new(-1.0, -1.999, 20, 20)
#new(-1.0, -1.000, 21, 21)
#new(-1.0, -1.499, 21, 21)
#new(-1.0, -1.999, 21, 21)


new_test(1000.0, 1000.000, 21, 21)
new_test(1000.0, 1000.499, 21, 21)
new_test(1000.0, 1000.999, 21, 21)
new_test(-1.0, -1.000, 20, 20)
new_test(-1.0, -1.499, 20, 20)
new_test(-1.0, -1.999, 20, 20)
new_test(-1.0, -1.000, 21, 21)
new_test(-1.0, -1.222, 21, 21)
new_test(-1.0, -1.499, 21, 21)
new_test(-1.0, -1.999, 21, 21)
