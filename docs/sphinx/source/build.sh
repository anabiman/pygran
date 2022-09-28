if [ -f _build ]; then rm -r _build; fi
if [ -f build ]; then rm -r build; fi

make html

cd _build/html

for d in _* ; do
    nd=${d//_}
    grep -rl "$d" *.html | xargs sed -i "s/${d}/${nd}/g"
    mv $d $nd
done

cd ../..

grep -rl '_autosummary' build/html | xargs sed -i 's/_autosummary/autosummary/g'
grep -rl '2014 drove.io' build/html | xargs sed -i 's/2014 drove.io/2022 PyGran/g'
grep -rl 'connectical' build/html | xargs sed -i 's/connectical/DEM open-source/g'

cp -r _build/html/* .
