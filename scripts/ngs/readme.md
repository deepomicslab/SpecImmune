
## Instructions

### Step 1: Create the Environment

Before using SpecImmune, ensure the required environment is set up. Please refer to the [SpecHLA environment setup guide](https://github.com/deepomicslab/SpecHLA) for detailed instructions.


move the bin directory to `ngs/bin` then run

`export LD_LIBRARY_PATH=path_to_your_spechla/spechla_env/lib/`

Note:
- SpecImmune-ngs requires a License of Novoalign in bin/. The License file of Novoalign should be put in the bin/ folder before next step.

---

### Step 2: Build the Database
The database files can be obtained from (https://drive.google.com/file/d/1mjivN6CllvO5myTorbQQv8zb-6Bop8tx/view?usp=sharing)

Build the database required by SpecImmune using the `script/makedb.sh` script.

Run the following command:

```bash
bash script/makedb.sh
```

### Step 3: Run the Script for Focus Genes
Run the script to analyze the target genes. Each gene has specific usage instructions, which can be found in the example.sh file.

Run the following command as an example:

```bash
bash example.sh
```