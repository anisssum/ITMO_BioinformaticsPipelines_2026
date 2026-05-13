process plot_coverage {
    publishDir "${params.outdir}/coverage", mode: 'copy'

    input:
    path bam  // ✅ только BAM файл

    output:
    path "coverage.png"
    path "depth.txt"

    script:
    """
    # Извлекаем sample_id из имени BAM файла
    SAMPLE_ID=\$(basename ${bam} .bam)
    
    samtools depth ${bam} > depth.txt
    
    python3 -c "
import matplotlib.pyplot as plt
with open('depth.txt') as f:
    depths = [int(l.split()[2]) for l in f]
plt.figure(figsize=(12,4))
plt.plot(depths, color='darkred', linewidth=0.5)
plt.title(f'Coverage Depth: \${SAMPLE_ID}')
plt.xlabel('Position')
plt.ylabel('Depth')
plt.grid(True, alpha=0.6)
plt.tight_layout()
plt.savefig('coverage.png', dpi=300)
    "
    """
}