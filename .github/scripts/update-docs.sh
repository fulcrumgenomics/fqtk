#!/bin/bash

set -euo pipefail

echo -e "<!-- start usage -->\n\`\`\`\`console\n" > usage.txt;
./target/debug/fqtk demux --help | sed -e 's_^[ ]*$__g' >> usage.txt
echo -e "\`\`\`\`\n<!-- end usage -->" >> usage.txt;
sed -e '/<!-- start usage -->/,/<!-- end usage -->/!b' -e '/<!-- end usage -->/!d;r usage.txt' -e 'd' README.md > README.md.new;
mv README.md.new README.md;
rm usage.txt;
