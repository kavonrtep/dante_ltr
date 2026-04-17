#!/bin/sh
# Install dante_ltr into $PREFIX. noarch package; just copies the source
# tree into $PREFIX/share/dante_ltr and symlinks executables into $PREFIX/bin.
set -x -e

DANTE_LTR_DIR="${PREFIX}/share/dante_ltr"

mkdir -p "${PREFIX}/bin"
mkdir -p "${DANTE_LTR_DIR}"
cp -r . "${DANTE_LTR_DIR}"

# Executables
ln -s "${DANTE_LTR_DIR}/dante_ltr"                 "${PREFIX}/bin/dante_ltr"
ln -s "${DANTE_LTR_DIR}/dante_ltr_solo"            "${PREFIX}/bin/dante_ltr_solo"
ln -s "${DANTE_LTR_DIR}/dante_ltr_to_library"      "${PREFIX}/bin/dante_ltr_to_library"
ln -s "${DANTE_LTR_DIR}/dante_ltr_summary"         "${PREFIX}/bin/dante_ltr_summary"
ln -s "${DANTE_LTR_DIR}/dante_reclassify"          "${PREFIX}/bin/dante_reclassify"
ln -s "${DANTE_LTR_DIR}/dante_ltr_gff3_to_canonical" "${PREFIX}/bin/dante_ltr_gff3_to_canonical"
ln -s "${DANTE_LTR_DIR}/clean_ltr.R"               "${PREFIX}/bin/clean_ltr.R"
