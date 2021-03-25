cd src || exit 1

if [[ $1 == test ]]
then
  make tests && ./tests
else
  make && ./main
fi