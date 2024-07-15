#!/bin/bash

# Create activate.d and deactivate.d directories
mkdir -p $CONDA_PREFIX/envs/SAMSA2/etc/conda/activate.d
mkdir -p $CONDA_PREFIX/envs/SAMSA2/etc/conda/deactivate.d

# Create activation script
cat > $CONDA_PREFIX/envs/SAMSA2/etc/conda/activate.d/env_vars.sh << EOL
#!/bin/bash
export ORIGINAL_PATH=\$PATH
export PATH=\${PATH}:$(dirname "$(realpath "$0")")
EOL

chmod +x ${CONDA_PREFIX}/envs/SAMSA2/etc/conda/activate.d/env_vars.sh

# Create deactivation script
cat > ${CONDA_PREFIX}/envs/SAMSA2/etc/conda/deactivate.d/unset_env_vars.sh << EOL
#!/bin/bash
export PATH=\$ORIGINAL_PATH
unset ORIGINAL_PATH
EOL

chmod +x ${CONDA_PREFIX}/envs/SAMSA2/etc/conda/deactivate.d/unset_env_vars.sh
