# network_generator
Creates a network of mouse gene interactions represented by their orthologous human genes interactions

## Set-up & Requirements
1. Install requirements from `requirements.txt`
2. If `./data/publications.csv` is absent, run `prep_work.py` to generate it. This is a one-time call to format all publication counts for each human gene.

## Demonstration
1. Load the website at [Gene Network Generator Site](http://bfx3.aap.jhu.edu/mli186/project/network.html) or click the url: http://bfx3.aap.jhu.edu/mli186/project/network.html
<kbd>![alt text](https://github.com/Meganmli/network_generator/blob/main/img/step1.png)</kbd>

2. Submit a mouse gene name. Autocompletion provides gene name options. Once submitted, be patient as the network takes time to generate.
<kbd>![alt text](https://github.com/Meganmli/network_generator/blob/main/img/step2.png)</kbd>

3. Explore the network by dragging around the nodes. Explore the static interaction table below.
<kbd>![alt text](https://github.com/Meganmli/network_generator/blob/main/img/step3.jpg)</kbd>