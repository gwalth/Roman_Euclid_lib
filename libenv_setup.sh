export ROMAN_EUCLID_LIB=$(cd `dirname $0` && pwd)
echo $ROMAN_EUCLID_LIB

export PATH=${ROMAN_EUCLID_LIB}/bin:${PATH}
export PYTHONPATH=${ROMAN_EUCLID_LIB}/lib:${PYTHONPATH}


