grep -rl '_static' build/html | xargs sed -i 's/_static/static/g'
grep -rl '_images' build/html | xargs sed -i 's/_images/images/g'
grep -rl '2014 drove.io' build/html | xargs sed -i 's/2014 drove.io/2019 PyGran/g'

mv build/html/_static build/html/static
mv build/html/_images build/html/images

