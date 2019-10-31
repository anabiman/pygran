grep -rl '_static' build/html | xargs sed -i 's/_static/static/g'
grep -rl '_images' build/html | xargs sed -i 's/_images/images/g'
grep -rl '_modules' build/html | xargs sed -i 's/_modules/modules/g'
grep -rl '_sources' build/html | xargs sed -i 's/_sources/sources/g'

grep -rl '_autosummary' build/html | xargs sed -i 's/_autosummary/autosummary/g'
grep -rl '2014 drove.io' build/html | xargs sed -i 's/2014 drove.io/2019 PyGran/g'

mv build/html/_static build/html/static
mv build/html/_images build/html/images
mv build/html/_autosummary build/html/autosummary
mv build/html/_sources build/html/sources
mv build/html/_modules build/html/modules
