if [ -z "$1" ]; then
    echo "Usage: $0 file.tex"
    exit 1
fi

FILE="$1"
DIR=$(dirname "$FILE")
BASE=$(basename "$FILE" .tex)

echo "Compiling $FILE..."

pdflatex -interaction=nonstopmode -output-directory="$DIR" "$FILE" > /dev/null

echo "Done. Output: $DIR/$BASE.pdf"
