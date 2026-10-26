#!/usr/bin/env python3
# 原第1行说明：指定在Unix环境中使用PATH中的python3解释器。通过python文件名启动时这行不决定conda环境。
# 原第2行：模块说明：这是连接官方RAISING的数据适配脚本，用于探索性扫描，不是已校准的选择检验。这里Official RAISING adapter不是说本文件是官方代码。
"""Official RAISING adapter. Exploratory scan, not calibrated selection inference."""
# 原第3行：导入命令行参数解析工具，读取prepare/tune/scan等子命令和参数。
import argparse
# 原第4行：导入上下文管理工具；后面使用redirect_stdout将训练的标准输出写入日志。
import contextlib
# 原第5行：导入CSV/TSV写入工具。
import csv
# 原第6行：导入gzip压缩文件读取工具。
import gzip
# 原第7行：导入哈希工具，用来给输入文件生成SHA256指纹。
import hashlib
# 原第8行：导入JSON读写工具。
import json
# 原第9行：导入操作系统接口，处理工作目录、文件大小、文件替换等。
import os
# 原第10行：导入Path，以对象方式拼接和操作文件路径。
from pathlib import Path
# 原第11行：导入sys；这个交付版本实际上没有用到它，可以保留但并非必需。
import sys
# 原第12行：导入NumPy，简称np，处理数值数组、缺失值、矩阵与随机数。
import numpy as np
# 原第13行：导入pandas，简称pd，处理样本表、位点表与模型目标表。
import pandas as pd

# 原第15行：定义辅助函数：把value保存成JSON文件。
def save_json(path, value):
    # 原第16行：先将value序列化为JSON文本，缩进2格并保留中文；再写入指定文件。
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False))

# 原第18行：定义辅助函数：计算文件SHA256指纹。
def digest(path):
    # 原第19行：将文件内容读成字节，计算哈希并返回十六进制字符串。实际主要用于小型元信息文件，不是给300万SNP大矩阵反复计算哈希。
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

# 原第21行：定义统一的文本文件打开函数。
def op(path):
    # 原第22行：文件名以.gz结尾时使用gzip文本模式读取，否则正常打开；不支持BCF。
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)

# 原第24行：定义读取VCF样本名的函数。
def header(path):
    # 原第25行：打开VCF文件，退出with块时自动关闭。
    with op(path) as f:
        # 原第26行：逐行读取，不一次把VCF全部装进内存。
        for line in f:
            # 原第27行：查找以#CHROM开头的VCF列标题行。
            if line.startswith('#CHROM'):
                # 原第28行：移除换行并按制表符拆分；Python从0计数，前9列是固定字段，第10列起是样本。
                names = line.rstrip('\n\r').split('\t')[9:]
                # 原第29行：比较样本名列表长度和去重后的长度，检查重复。
                if len(names) != len(set(names)):
                    # 原第30行：存在重复样本名时停止，避免一个ID对应多列。
                    raise ValueError('Duplicate VCF sample IDs')
                # 原第31行：找到合法标题后返回样本名列表。
                return names
    # 原第32行：读完整个文件仍找不到标题时停止。
    raise ValueError('No #CHROM header')

# 原第34行：定义逐位点生成器；wanted是希望保留且按指定顺序排列的样本名，max_records=0表示不限制读取记录数。
def records(path, wanted, max_records=0):
    # 原第35行：读取VCF中的全部样本名。
    names = header(path)
    # 原第36行：找出wanted中没有出现在VCF里的样本。
    missing = set(wanted) - set(names)
    # 原第37行：如果缺少任何希望读取的样本，进入报错分支。
    if missing:
        # 原第38行：列出缺失样本并停止，避免静默改变样本集合。
        raise ValueError(f'VCF lacks samples: {sorted(missing)}')
    # 原第39行：找到每个保留样本在VCF数据行中的列索引；加9是跳过固定列。
    indices = [names.index(s) + 9 for s in wanted]
    # 原第40行：打开VCF，开始逐行扫描。
    with op(path) as f:
        # 原第41行：初始化已遇到的数据记录数。
        n = 0
        # 原第42行：逐行遍历VCF。
        for line in f:
            # 原第43行：判断是否为元信息或标题行。
            if line.startswith('#'):
                # 原第44行：跳过所有以#开头的行。
                continue
            # 原第45行：当前数据记录计数加1，注意统计的是过滤前记录数。
            n += 1
            # 原第46行：如果设置了max_records且已超过限制，停止读取。
            if max_records and n > max_records:
                # 原第47行：退出读取循环。smoke最多检查前20000条数据记录，而不是保证保留20000个SNP。
                break
            # 原第48行：将当前位点按制表符拆分成字段。
            r = line.rstrip('\n\r').split('\t')
            # 原第49行：检查行宽是否等于9个固定字段加样本数。
            if len(r) != 9 + len(names):
                # 原第50行：行字段数异常时停止。
                raise ValueError(f'Invalid VCF width, record {n}')
            # 原第51行：只接受REF和ALT均为单个A/C/G/T的双等位SNP，以及FILTER为PASS或.的记录。FILTER=.仅表示没有过滤信息，不等于通过全部质量检查。
            if len(r[3]) != 1 or r[3] not in 'ACGT' or len(r[4]) != 1 or r[4] not in 'ACGT' or r[6] not in ('.','PASS'):
                # 原第52行：跳过不符合这些条件的位点；因此当前函数不适用于直接扫描SV或多等位位点。
                continue
            # 原第53行：拆分FORMAT字段，例如GT:DP:GQ。
            fmt = r[8].split(':')
            # 原第54行：检查FORMAT是否包含GT基因型。
            if 'GT' not in fmt:
                # 原第55行：没有GT则停止。
                raise ValueError('GT absent')
            # 原第56行：记录GT在每个样本冒号分隔字段中的位置。
            gi = fmt.index('GT')
            # 原第57行：为所有保留样本分配float32数组，先全部设为NaN；这是当前一个位点的剂量向量。
            x = np.full(len(wanted), np.nan, dtype=np.float32)
            # 原第58行：遍历保留样本；j是输出数组中的位置，i是VCF行里的列位置。
            for j, i in enumerate(indices):
                # 原第59行：拆分该样本字段。
                fields = r[i].split(':')
                # 原第60行：取出GT；字段不足则作为缺失值。
                g = fields[gi] if gi < len(fields) else '.'
                # 原第61行：将0|1等分相写法转换成0/1，再拆成两个等位基因；这里不保留相位信息。
                alleles = g.replace('|','/').split('/')
                # 原第62行：任一等位基因缺失时，整个样本在该位点都按缺失处理。
                if '.' in alleles:
                    # 原第63行：跳过该样本，本位点剂量仍保持NaN。
                    continue
                # 原第64行：要求恰好两个等位基因且只能为0或1；当前流程假设二倍体。
                if len(alleles) != 2 or any(v not in ('0','1') for v in alleles):
                    # 原第65行：发现单倍体、多倍体或非法编码时停止。
                    raise ValueError(f'Non-diploid/non-biallelic GT: {g}')
                # 原第66行：将GT转换为ALT剂量：0/0→0，0/1或1/0→1，1/1→2。
                x[j] = sum(int(v) for v in alleles)
            # 原第67行：计算该位点在保留样本中的缺失比例。
            miss = float(np.isnan(x).mean())
            # 原第68行：缺失率超过10%或完全缺失时不保留。
            if miss > 0.1 or miss == 1:
                # 原第69行：跳过该位点。
                continue
            # 原第70行：非缺失剂量均值除以2，得到ALT等位基因频率。
            af = float(np.nanmean(x) / 2)
            # 原第71行：较小的等位基因频率是min(AF,1-AF)。
            maf = min(af, 1-af)
            # 原第72行：判断MAF是否小于0.05；这是匹配样本集合内重新计算的MAF。
            if maf < 0.05:
                # 原第73行：跳过低频位点。
                continue
            # 原第74行：用2×AF，也就是观测剂量均值，填补剩余缺失；可能得到0.6等非整数剂量。
            x[np.isnan(x)] = 2*af
            # 原第75行：逐条产出原字段、剂量、ALT频率、MAF、填补前缺失率。yield让外层逐位点处理，避免整份VCF同时载入。
            yield r, x, af, maf, miss

# 原第77行：定义prepare子命令的处理函数；a是命令行参数对象。
def prepare(a):
    # 原第78行：把输出目录转换为Path对象。
    out = Path(a.out)
    # 原第79行：检查输出目录是否已存在。
    if out.exists():
        # 原第80行：已存在则停止，避免覆盖旧的预处理结果；失败后留下目录也需要使用新目录或核查后清理。
        raise ValueError(f'{out} already exists; use a new directory to preserve provenance')
    # 原第81行：读取环境CSV。
    env = pd.read_csv(a.env)
    # 原第82行：将第一列命名为sample_id；所以要求原CSV第一列确实为材料ID。
    env = env.rename(columns={env.columns[0]:'sample_id'})
    # 原第83行：样本ID转成字符串并清理首尾空格。
    env['sample_id'] = env.sample_id.astype(str).str.strip()
    # 原第84行：检测环境表中的重复样本ID。
    if env.sample_id.duplicated().any():
        # 原第85行：重复则停止，不随意挑一行。
        raise ValueError('Duplicate environment sample IDs')
    # 原第86行：生成BIO1到BIO19的列名。range的右端20不包含在内。
    bio = [f'BIO{i}' for i in range(1,20)]
    # 原第87行：要求存在经纬度及全部19个BIO变量。
    required = ['lon','lat'] + bio
    # 原第88行：将这些列转换成数值，遇到无法转换的文本时直接报错。
    env[required] = env[required].apply(pd.to_numeric, errors='raise')
    # 原第89行：good标记环境数值均有限且经纬度在合法范围内的样本；不代表坐标来源或生态意义已经核实。
    good = np.isfinite(env[required]).all(axis=1) & env.lon.between(-180,180) & env.lat.between(-90,90)
    # 原第90行：读取VCF样本顺序。
    samples = header(a.vcf)
    # 原第91行：构建环境信息完整有效的样本ID集合。
    valid = set(env.loc[good,'sample_id'])
    # 原第92行：按VCF顺序保留同时有有效环境信息的样本，保证基因型与环境逐行对应。
    kept = [s for s in samples if s in valid]
    # 原第93行：检查匹配样本是否少于15份。
    if len(kept) < 15:
        # 原第94行：不足15份时停止；15是程序防护阈值，并不是神经网络统计功效充足的标准。
        raise ValueError('Fewer than 15 matched complete samples')
    # 原第95行：以ID为索引，按照kept顺序重排环境表，再恢复普通索引。
    matched = env.set_index('sample_id').loc[kept].reset_index()
    # 原第96行：创建输出目录及必要父目录。
    out.mkdir(parents=True)
    # 原第97行：保存保留样本的完整环境与坐标信息。
    matched.to_csv(out/'samples.csv', index=False)
    # 原第98行：初始化预测目标表y，先放入样本ID。
    y = pd.DataFrame({'sample_id':kept})
    # 原第99行：导入标准化工具，将不同量纲的环境变量变成均值0、标准差1。
    from sklearn.preprocessing import StandardScaler
    # 原第100行：导入主成分分析工具。
    from sklearn.decomposition import PCA
    # 原第101行：分别处理两个变量组：BIO1–BIO11为温度组；BIO12–BIO19为降水组。这是本适配脚本设置，不是官方RAISING强制要求。
    for label, columns in [('Temperature_PC1',bio[:11]),('Precipitation_PC1',bio[11:])]:
        # 原第102行：为当前组新建标准化器。
        scaler = StandardScaler()
        # 原第103行：用全部匹配样本估计均值/尺度并标准化；正式外层预测验证时应改为训练折拟合。
        z = scaler.fit_transform(matched[columns])
        # 原第104行：只提取一个主成分。
        pca = PCA(n_components=1)
        # 原第105行：拟合PCA，并取所有样本在PC1上的得分。
        scores = pca.fit_transform(z)[:,0]
        # 原第106行：取PC1对各环境变量的载荷。
        loading = pca.components_[0]
        # 原第107行：令绝对值最大载荷为正，固定PC符号；这不保证PC越大就是越热/越湿。
        sign = 1 if loading[np.argmax(abs(loading))] >= 0 else -1
        # 原第108行：同时翻转样本得分的符号，保存到对应目标列。
        y[label] = scores*sign
        # 原第109行：保存当前组的变量名、标准化均值/尺度及符号调整后的载荷。
        pd.DataFrame({'variable':columns,'mean':scaler.mean_,'scale':scaler.scale_,'loading':loading*sign}).to_csv(out/f'{label}_loadings.csv',index=False)
        # 原第110行：保存PC1解释的方差比例；这是对应变量组的方差，不是网络预测准确率。
        save_json(out/f'{label}_variance.json',{'explained_variance_ratio':float(pca.explained_variance_ratio_[0])})
    # 原第111行：保存包含sample_id及两个目标列的targets.csv。
    y.to_csv(out/'targets.csv',index=False)
    # 原第112行：初始化保留位点数。
    count = 0
    # 原第113行：用集合记录已遇到的染色体。
    seen_chrom = set()
    # 原第114行：初始化上一条记录的染色体。
    previous_chrom = None
    # 原第115行：初始化上一条记录的位置。
    previous_pos = -1
    # 原第116行：初始化上一条记录的变异标识。
    previous_key = None
    # 原第117行：同时打开二进制剂量矩阵文件和位点TSV。二进制存每个SNP的样本剂量，TSV存坐标和注释字段。
    with open(out/'dosage.float32','wb') as binary, open(out/'variants.tsv','w') as table:
        # 原第118行：建立制表符分隔的写入器。
        writer = csv.writer(table,delimiter='\t')
        # 原第119行：写入位点表的列名。index是过滤后全体SNP的0-based行号，pos仍是VCF的1-based基因组坐标。
        writer.writerow(['index','chrom','pos','variant_id','original_id','ref','alt','alt_frequency','maf','missing_rate'])
        # 原第120行：逐条取得质控后的保留位点。
        for r,x,af,maf,miss in records(a.vcf,kept,a.max_records):
            # 原第121行：取得染色体与整数位置。
            chrom,pos = r[0],int(r[1])
            # 原第122行：判断是否进入新染色体。
            if chrom != previous_chrom:
                # 原第123行：若当前染色体以前已出现，说明染色体记录被拆散。
                if chrom in seen_chrom:
                    # 原第124行：这种不连续染色体排列会破坏滑窗，直接停止。
                    raise ValueError('VCF chromosome blocks are not contiguous')
                # 原第125行：记录新染色体。
                seen_chrom.add(chrom)
                # 原第126行：重置染色体内的位置检查。
                previous_pos = -1
            # 原第127行：检查当前坐标是否小于上一条坐标。
            if pos < previous_pos:
                # 原第128行：位置逆序时停止。
                raise ValueError('VCF positions are not sorted')
            # 原第129行：构造chr:pos:ref:alt形式的完整变异标识，不完全依赖原VCF ID。
            key = f'{chrom}:{pos}:{r[3]}:{r[4]}'
            # 原第130行：检查是否与紧挨着的上一条变异完全重复。
            if key == previous_key:
                # 原第131行：相邻重复则停止；这里只检查相邻重复，并非全局ID重复扫描。
                raise ValueError(f'Duplicate adjacent variant: {key}')
            # 原第132行：将当前位点的样本剂量以float32连续写入二进制文件。
            x.tofile(binary)
            # 原第133行：将位点索引、坐标、变异ID、原ID、等位基因及QC信息写入TSV。
            writer.writerow([count,chrom,pos,key,r[2],r[3],r[4],af,maf,miss])
            # 原第134行：更新上一条记录的信息。
            previous_chrom,previous_pos,previous_key = chrom,pos,key
            # 原第135行：保留位点计数加1。
            count += 1
    # 原第136行：检查是否一个位点都没有保留。
    if count == 0:
        # 原第137行：全部被过滤时停止。
        raise ValueError('No variants passed filtering')
    # 原第138行：开始写入预处理元信息清单manifest.json。
    save_json(out/'manifest.json',{
        # 原第139行：记录样本数、SNP数以及二进制布局：行是SNP、列是样本。
        'samples':len(kept),'variants':count,'layout':'variants x samples, float32 dosage',
        # 原第140行：记录VCF绝对路径和字节大小。
        'vcf':str(Path(a.vcf).resolve()),'vcf_size':os.path.getsize(a.vcf),
        # 原第141行：记录VCF修改时间和环境CSV哈希；注意这里没有给整个VCF内容计算哈希。
        'vcf_mtime_ns':os.stat(a.vcf).st_mtime_ns,'environment_sha256':digest(a.env),
        # 原第142行：记录只有基因型、没有环境记录的样本。
        'vcf_only_samples':sorted(set(samples)-set(env.sample_id)),
        # 原第143行：记录只有环境、没有基因型的样本。
        'environment_only_samples':sorted(set(env.sample_id)-set(samples)),
        # 原第144行：记录环境缺失或坐标非法的材料。
        'invalid_environment_samples':env.loc[~good,'sample_id'].tolist(),
        # 原第145行：记录输入限制与MAF、缺失率阈值。
        'max_input_records':a.max_records,'maf_min':0.05,'missing_max':0.1,
        # 原第146行：明确写出全样本PCA与均值填补的限制：当前是探索性扫描，不是独立留出预测验证。
        'warning':'Environmental PCA and mean imputation use all matched samples. Exploratory scan only; not held-out predictive validation.'})
    # 原第147行：打印最终匹配样本和保留SNP数。
    print(f'Prepared {len(kept)} samples, {count} SNPs')

# 原第149行：定义读取准备数据的函数。
def load(data):
    # 原第150行：取得准备目录绝对路径，避免后面切换工作目录后找不到文件。
    d = Path(data).resolve()
    # 原第151行：读取manifest，了解矩阵形状。
    m = json.loads((d/'manifest.json').read_text())
    # 原第152行：把磁盘剂量文件映射为只读float32数组，形状为SNP数×样本数；实际按需访问，不一次性复制整矩阵。
    x = np.memmap(d/'dosage.float32',dtype=np.float32,mode='r',shape=(m['variants'],m['samples']))
    # 原第153行：读取位点表，染色体强制按字符串处理。
    v = pd.read_csv(d/'variants.tsv',sep='\t',dtype={'chrom':str})
    # 原第154行：读取目标表，删除样本ID，只留下两个数值目标列；依赖prepare已保证顺序一致。
    y = pd.read_csv(d/'targets.csv').drop(columns='sample_id')
    # 原第155行：返回目录、剂量矩阵、位点表及目标矩阵。
    return d,x,v,y

# 原第157行：定义滑窗生成器；width和step单位都是SNP个数，不是bp。
def windows(v,width,step):
    # 原第158行：检查步长合法性。
    if step < 1 or step > width:
        # 原第159行：要求步长为正且不大于窗口大小。
        raise ValueError('Require 1 <= step <= window')
    # 原第160行：按染色体分组，保留输入染色体顺序，避免窗口跨染色体。
    for chrom,g in v.groupby('chrom',sort=False):
        # 原第161行：取得该染色体在全局SNP表中的起止索引，右端为不包含的终点。
        start,end = int(g.index.min()),int(g.index.max())+1
        # 原第162行：通常用width；整条染色体SNP不足width时用该染色体实际长度。
        size = min(width,end-start)
        # 原第163行：按照step生成能放下完整窗口的起点。
        starts = list(range(start,end-size+1,step))
        # 原第164行：检查最后一个窗口是否到达染色体末端。
        if starts[-1] != end-size:
            # 原第165行：如果没有，补一个向前对齐的尾窗口；末端可能因此额外重叠。
            starts.append(end-size)
        # 原第166行：遍历窗口起点。
        for s in starts:
            # 原第167行：产出染色体、全局起点、右端不包含的终点；如0:1000包含索引0到999。
            yield chrom,s,s+size

# 原第169行：定义structure子命令，计算背景遗传PC。
def structure(a):
    # 原第170行：载入准备数据；这里有些返回值不会使用。
    d,x,v,y = load(a.data)
    # 原第171行：读取与准备数据一致的样本顺序。
    samples = pd.read_csv(d/'samples.csv').sample_id.tolist()
    # 原第172行：初始化样本×样本的关系矩阵累加器，125份样本就是125×125。
    gram = np.zeros((len(samples),len(samples)))
    # 原第173行：初始化背景SNP计数。
    count = 0
    # 原第174行：从传入的背景VCF逐位点读取；run_pepo.sh传入的是LD过滤集，但本函数自己不执行LD过滤。
    for r,g,af,maf,miss in records(a.vcf,samples):
        # 原第175行：标准化剂量：(g-2p)/sqrt(2p(1-p))。
        z = (g-2*af)/np.sqrt(2*af*(1-af))
        # 原第176行：计算当前位点标准化向量的外积，累加到关系矩阵。
        gram += np.outer(z,z)
        # 原第177行：背景SNP计数加1。
        count += 1
    # 原第178行：检查有效背景位点是否不足10个。
    if count < 10:
        # 原第179行：位点过少则停止。
        raise ValueError('Too few LD-pruned background SNPs')
    # 原第180行：除以背景位点数，得到平均关系矩阵。
    gram /= count
    # 原第181行：对实对称矩阵进行特征分解。
    eigen,vec = np.linalg.eigh(gram)
    # 原第182行：把特征值从大到小排序，选前10个。
    order = np.argsort(eigen)[::-1][:10]
    # 原第183行：用对应特征向量乘特征值平方根，得到缩放后的遗传轴；负的小数值误差截为0。该缩放不改变相关/线性回归的解释结果。
    pc = pd.DataFrame(vec[:,order]*np.sqrt(np.maximum(eigen[order],0)),columns=[f'genetic_PC{i+1}' for i in range(len(order))])
    # 原第184行：在PC表第一列加入sample_id。
    pc.insert(0,'sample_id',samples)
    # 原第185行：保存背景遗传PC。
    pc.to_csv(d/'background_genetic_PCs.csv',index=False)
    # 原第186行：记录背景SNP数与文件路径，明确这些PC只作诊断，没有加入官方网络做校正。
    save_json(d/'background_structure.json',{'snp_count':count,'vcf':str(Path(a.vcf).resolve()),'warning':'Diagnostics only: official RAISING does not condition on these PCs in this adapter.'})
    # 原第187行：打印背景位点数。
    print(f'Background structure: {count} SNPs')

# 原第189行：定义TensorFlow随机初始化配置函数。
def setup_tf(seed):
    # 原第190行：仅需要网络时导入TensorFlow，prepare等步骤不需要它。
    import tensorflow as tf
    # 原第191行：设置Python/NumPy/TensorFlow相关随机种子；官方feature_importance内部还会重新设置其自己的种子。
    tf.keras.utils.set_random_seed(seed)
    # 原第192行：返回TensorFlow模块。
    return tf

# 原第194行：定义调参架构函数。a是参数，d是准备目录，x/v/y是剂量、位点、目标，width是网络输入SNP数。
def architecture(a,d,x,v,y,width):
    # 原第195行：导入TensorFlow。
    import tensorflow as tf
    # 原第196行：这里才导入RAISING官方的hp_optimization调参函数；之前的预处理是自定义代码。
    from RAISING.hp_tune_select import hp_optimization
    # 原第197行：取得结果目录绝对路径。
    root = Path(a.out).resolve()
    # 原第198行：为该输入长度创建调参目录，例如tune_w1000。
    dest = root/f'tune_w{width}'
    # 原第199行：创建目录，已存在则沿用，后面再核对参数。
    dest.mkdir(parents=True,exist_ok=True)
    # 原第200行：指定官方函数保存模型的位置。文件包含模型结构与训练后的权重；后续DeepFeatImp调用会重置权重再训练。
    model = dest/'architecture.keras'
    # 原第201行：构造调参元信息，记录窗口长度、数据与目标指纹、随机种子、尝试数和训练上限。
    meta = {'width':width,'data_manifest_sha256':digest(d/'manifest.json'),'targets_sha256':digest(d/'targets.csv'),'seed':a.seed,'trials':a.trials,'epochs':a.epochs}
    # 原第202行：若模型已经存在，尝试复用。
    if model.exists():
        # 原第203行：读取原run.json并比较全部元信息。
        if json.loads((dest/'run.json').read_text()) != meta:
            # 原第204行：数据或设置不一致时停止，不把旧模型混用到新任务。
            raise ValueError('Tuning configuration changed; choose a fresh --out directory')
        # 原第205行：一致则返回现有模型路径，不重新调参。
        return model
    # 原第206行：枚举能容纳当前长度的调参候选窗口；这里step=width，基本不重叠，尾部仍可能补窗。
    candidates = [(s,e) for c,s,e in windows(v,width,width) if e-s==width]
    # 原第207行：检查是否存在合法候选窗口。
    if not candidates:
        # 原第208行：没有则停止。
        raise ValueError(f'No chromosome has {width} SNPs')
    # 原第209行：用固定随机种子随机选择一个候选窗口；当前调参只用这一个窗口。
    s,e = candidates[np.random.default_rng(a.seed).integers(len(candidates))]
    # 原第210行：开始定义自定义超参数搜索空间，不是复制论文全部参数。
    config = {
        # 原第211行：隐藏层数1–2，每层神经元数8–32；当前官方构造方式同一次模型各隐藏层使用同一个units值。
        'layers':{'min_value':1,'max_value':2},'units':{'min_value':8,'max_value':32},
        # 原第212行：学习率在0.0001到0.001之间按对数尺度搜索。
        'learning_rate':{'min_value':0.0001,'max_value':0.001,'sampling':'log'},
        # 原第213行：激活函数固定ReLU，L1正则化强度在1e-5到1e-3之间搜索。
        'activation':{'values':['relu']},'l1_regularizer':{'min_value':0.00001,'max_value':0.001,'sampling':'log'},
        # 原第214行：L2正则化强度在1e-5到1e-3之间搜索。
        'l2_regularizer':{'min_value':0.00001,'max_value':0.001,'sampling':'log'},
        # 原第215行：dropout比例从0.1到0.5搜索。
        'dropout':{'min_value':0.1,'max_value':0.5,'sampling':'linear'},
        # 原第216行：优化器固定Adam，权重初始化固定he_normal，结束配置字典。
        'optimizer':{'values':['adam']},'weight_initializer':{'values':['he_normal']}}
    # 原第217行：保存搜索空间供复查。
    save_json(dest/'search_space.json',config)
    # 原第218行：保存原工作目录。
    old = Path.cwd()
    # 原第219行：进入try/finally，确保即使失败也恢复工作目录。
    try:
        # 原第220行：切换到本次调参目录，官方函数产生的相对路径文件会保存在这里。
        os.chdir(dest)
        # 原第221行：打开训练日志并重定向标准输出；没有设置实时刷新，所以日志可能暂时为空。stderr不在这里重定向。
        with open('training.log','w') as log, contextlib.redirect_stdout(log):
            # 原第222行：调用官方调参：取选中窗口x[s:e]并转置，输入从SNP×样本变为样本×SNP；y为样本×两个环境目标，目标类型连续，优化val_loss。
            hp_optimization(np.asarray(x[s:e].T),y,output_class='continuous',objective_fun='val_loss',
                # 原第223行：指定贝叶斯超参数优化、搜索空间JSON和模型保存路径；贝叶斯优化不等于BayPass的贝叶斯关联检验。
                algorithm='Bayesian',config_file=str(dest/'search_space.json'),model_file=str(model),
                # 原第224行：设置最多trials种配置、每次最多epochs轮，批量16，early stopping耐心10，最小改善1e-4。
                max_trials=a.trials,epochs=a.epochs,batch_size=16,patience=10,min_delta=0.0001,
                # 原第225行：随机留25%样本作内部验证；固定随机种子，verbose=0不逐轮打印。未按地点/亲缘分组；官方内部先标准化再划分，不能作为严格外部验证。
                test_size=0.25,random_state=a.seed,seed=a.seed,verbose=0)
        # 原第226行：官方调用成功后保存参数与数据指纹。
        save_json(dest/'run.json',meta)
        # 原第227行：保存调参窗口的全局起止索引。
        save_json(dest/'tuning_window.json',{'start_index':s,'end_index_exclusive':e})
    # 原第228行：不论成功或失败，都进入清理分支。
    finally:
        # 原第229行：恢复原工作目录。
        os.chdir(old)
        # 原第230行：清理Keras会话，释放模型相关资源。
        tf.keras.backend.clear_session()
    # 原第231行：返回模型文件路径。
    return model

# 原第233行：定义tune子命令入口。
def tune(a):
    # 原第234行：载入准备数据。
    d,x,v,y = load(a.data)
    # 原第235行：设置随机种子。
    setup_tf(a.seed)
    # 原第236行：从实际窗口中收集不同长度并去重排序；通常只有1000，若某染色体不足1000则为短长度单独调参。
    for width in sorted({e-s for c,s,e in windows(v,a.window,a.window)}):
        # 原第237行：打印当前调参长度，flush=True立即输出，便于看状态。
        print('Tuning width',width,flush=True)
        # 原第238行：为该长度调用上面的architecture函数。
        architecture(a,d,x,v,y,width)

# 原第240行：定义scan子命令入口。
def scan(a):
    # 原第241行：导入RAISING官方feature_importance函数，真实的重复训练及DeepFeatImp计算在官方函数内部。
    from RAISING.hp_tune_select import feature_importance
    # 原第242行：设置随机种子并取得TensorFlow模块。
    tf = setup_tf(a.seed)
    # 原第243行：载入准备数据。
    d,x,v,y = load(a.data)
    # 原第244行：取得输出目录绝对路径。
    root = Path(a.out).resolve()
    # 原第245行：创建目录。
    root.mkdir(parents=True,exist_ok=True)
    # 原第246行：开始记录扫描配置，首先记录数据和环境目标指纹。
    cfg = {'manifest_sha256':digest(d/'manifest.json'),'targets_sha256':digest(d/'targets.csv'),
           # 原第247行：记录窗口、步长、训练轮数、重复数和随机种子。
           'window':a.window,'step':a.step,'epochs':a.epochs,'repeats':a.repeats,'seed':a.seed}
    # 原第248行：定义扫描配置文件位置。
    configpath = root/'scan_config.json'
    # 原第249行：如果旧配置存在，比较本次是否一致。
    if configpath.exists() and json.loads(configpath.read_text()) != cfg:
        # 原第250行：不一致则停止，避免不同设置的窗口文件混在一起。
        raise ValueError('Scan configuration changed; use a fresh --out')
    # 原第251行：保存扫描配置。
    save_json(configpath,cfg)
    # 原第252行：初始化本次遍历窗口计数；它包括遇到后跳过的已完成窗口，并非本次新训练数量。
    completed = 0
    # 原第253行：按染色体依次生成滑窗。当前是串行窗口循环，不是多窗口并行。
    for chrom,s,e in windows(v,a.window,a.step):
        # 原第254行：若设置了limit_windows且已经遍历到上限，停止。
        if a.limit_windows and completed >= a.limit_windows:
            # 原第255行：退出扫描循环；limit=2用于试运行。
            break
        # 原第256行：本次遍历计数加1。
        completed += 1
        # 原第257行：按本窗口输入长度找到对应调参模型。
        model = root/f'tune_w{e-s}'/'architecture.keras'
        # 原第258行：检查模型是否存在。
        if not model.exists():
            # 原第259行：模型缺失则提示先运行tune。
            raise ValueError(f'Missing {model}. Run tune first.')
        # 原第260行：读取模型配套元信息。
        metadata = json.loads((model.parent/'run.json').read_text())
        # 原第261行：核对模型对应的数据与目标是否与本次一致。
        if metadata['data_manifest_sha256'] != digest(d/'manifest.json') or metadata['targets_sha256'] != digest(d/'targets.csv'):
            # 原第262行：不一致则停止。
            raise ValueError('Model and prepared data mismatch')
        # 原第263行：为当前窗口生成补零后的目录名，使用全局SNP索引，不是基因组bp。
        dest = root/'windows'/f'w_{s:09d}_{e:09d}'
        # 原第264行：创建本窗口目录。
        dest.mkdir(parents=True,exist_ok=True)
        # 原第265行：定义窗口最终重要性文件路径。
        final = dest/'importance.csv'
        # 原第266行：如果该文件已经存在，认为该窗口完成。
        if final.exists():
            # 原第267行：打印断点跳过信息。
            print('Resume: skipping',dest.name,flush=True)
            # 原第268行：跳到下个窗口；这里主要依赖文件存在，完整性还需后处理核查。
            continue
        # 原第269行：打印将要训练的染色体与索引范围。
        print('Training',chrom,s,e,flush=True)
        # 原第270行：保存当前目录。
        old = Path.cwd()
        # 原第271行：开始带finally的训练过程。
        try:
            # 原第272行：切到窗口目录，隔离官方函数的固定文件名输出。
            os.chdir(dest)
            # 原第273行：将标准输出写入该窗口training.log，同样可能存在缓冲。
            with open('training.log','w') as log, contextlib.redirect_stdout(log):
                # 原第274行：调用官方feature_importance，输入为当前窗口的样本×SNP剂量和两个环境目标。
                result = feature_importance(np.asarray(x[s:e].T),y,
                    # 原第275行：提供与每一列剂量一一对应的变异ID，以及调参得到的模型文件。
                    feature_set=v.iloc[s:e].variant_id.tolist(),model_file=str(model),
                    # 原第276行：iteration指定每窗重复训练次数，方法固定DeepFeatImp，目标为连续变量。当前全量设置iteration=5。
                    iteration=a.repeats,feature_method='DeepFeatImp',output_class='continuous',
                    # 原第277行：指定官方函数产生的平均权重模型文件名。
                    train_model_file=str(dest/'upstream_average_model.h5'),
                    # 原第278行：每次训练最多epochs轮、批量16、提前停止参数、关闭逐轮打印。DeepFeatImp分支在全体样本上训练并监控训练loss。
                    epochs=a.epochs,batch_size=16,patience=10,min_delta=0.0001,verbose=0,
                    # 原第279行：给官方内部划分参数提供随机种子；DeepFeatImp分支不使用独立测试集评估，不能据此声称泛化准确。
                    random_state=a.seed)
            # 原第280行：取出两个目标的重要性列并转换为数值。
            vals = result[y.columns].astype(float)
            # 原第281行：检查输出SNP数是否正确，分数是否都为有限数。
            if len(result) != e-s or not np.isfinite(vals.to_numpy()).all():
                # 原第282行：输出异常则停止，不将非法结果视为成功。
                raise ValueError('Invalid official feature importance output')
            # 原第283行：先写临时结果文件。
            result.to_csv(dest/'importance.tmp',index=False)
            # 原第284行：完成后用os.replace将临时文件替换成最终文件名，降低半写入文件被认作完成的风险。
            os.replace(dest/'importance.tmp',final)
            # 原第285行：说明官方实现会对多次初始化训练的权重做平均。
            # Upstream averages weights across independently initialized networks.
            # 原第286行：我们保留其逐次训练后平均的重要性，不把平均权重模型当成可靠预测器。
            # Retain importance only; do not use averaged weights for prediction.
            # 原第287行：删除不用于预测的官方平均权重模型，节省磁盘。官方产生的其他固定CSV/日志仍可能留在窗口目录。
            (dest/'upstream_average_model.h5').unlink(missing_ok=True)
            # 原第288行：保存窗口染色体、索引和架构文件哈希；若在写importance.csv后、写此文件前中断，需修复/重跑该窗口。
            save_json(dest/'window.json',{'chrom':chrom,'start':s,'end':e,'architecture_sha256':digest(model)})
        # 原第289行：不论训练成功与否，执行清理。
        finally:
            # 原第290行：恢复原工作目录。
            os.chdir(old)
            # 原第291行：清理Keras会话。
            tf.keras.backend.clear_session()

# 原第293行：定义aggregate子命令，汇总窗口重要性；它没有计算尾部概率或Matrix Profile。
def aggregate(a):
    # 原第294行：载入数据；这里主要使用v和y。
    d,x,v,y = load(a.data)
    # 原第295行：取得输出目录。
    root = Path(a.out)
    # 原第296行：分配SNP数×目标数矩阵，累加每个SNP的重要性。
    total = np.zeros((len(v),len(y.columns)))
    # 原第297行：分配同形状矩阵，累加重要性的平方，用于计算跨窗波动。
    total2 = np.zeros_like(total)
    # 原第298行：初始化每个位点的窗口覆盖次数。
    counts = np.zeros(len(v),dtype=int)
    # 原第299行：读取扫描配置。
    config = json.loads((root/'scan_config.json').read_text())
    # 原第300行：检查准备数据指纹是否与扫描时一致。
    if config['manifest_sha256'] != digest(d/'manifest.json'):
        # 原第301行：不一致则停止。
        raise ValueError('Prepared data changed')
    # 原第302行：查找所有窗口的importance.csv并按路径排序。
    for f in sorted((root/'windows').glob('*/importance.csv')):
        # 原第303行：读取当前窗口重要性。
        z = pd.read_csv(f)
        # 原第304行：读取该窗口的索引元信息。
        m = json.loads((f.parent/'window.json').read_text())
        # 原第305行：取得全局起止索引。
        s,e = m['start'],m['end']
        # 原第306行：严格比较窗口变异ID顺序与全局位点表切片顺序。
        if z.features.tolist() != v.iloc[s:e].variant_id.tolist():
            # 原第307行：顺序不一致则停止，防止把分数赋给错误SNP。
            raise ValueError('Importance/variant order mismatch')
        # 原第308行：将两个环境目标的重要性列变成数值数组。
        arr = z[y.columns].to_numpy(dtype=float)
        # 原第309行：注释表达的是尺寸缩放的意图；并不保证不同窗口的重要性因此完全可比。
        # Size normalization makes small-chromosome models less incomparable.
        # 原第310行：将官方重要性乘窗口SNP数。若1000个位点的重要性近似均分且总和1，平均分由0.001变到约1；这不是P值转换。
        arr *= e-s
        # 原第311行：将本窗口缩放后的分数累加到对应全局位点。
        total[s:e] += arr
        # 原第312行：同时累加分数平方。
        total2[s:e] += arr**2
        # 原第313行：这些位点的覆盖窗口数加1。
        counts[s:e] += 1
    # 原第314行：只保留至少有一个窗口结果的位点，因此部分扫描也能生成汇总表。
    keep = counts > 0
    # 原第315行：复制这些位点的坐标与QC信息。
    result = v.loc[keep].copy()
    # 原第316行：加入窗口覆盖次数n_windows。
    result['n_windows'] = counts[keep]
    # 原第317行：逐个环境目标汇总。
    for j,c in enumerate(y.columns):
        # 原第318行：累加值除覆盖次数，得到跨窗口平均重要性。
        mean = total[keep,j]/counts[keep]
        # 原第319行：按sqrt(E[x²]-E[x]²)计算跨窗总体标准差；截掉浮点负误差。只有一个窗口时为0，不代表稳定性已验证。
        sd = np.sqrt(np.maximum(0,total2[keep,j]/counts[keep]-mean**2))
        # 原第320行：写入该目标的平均重要性列。
        result[c+'_importance_scaled'] = mean
        # 原第321行：写入跨窗口标准差列。
        result[c+'_window_sd'] = sd
        # 原第322行：按重要性从高到低排名；并列取最小名次。排名1最高，不是显著性检验。
        result[c+'_rank'] = pd.Series(mean).rank(ascending=False,method='min').to_numpy()
    # 原第323行：保存压缩TSV全局重要性表。
    result.to_csv(root/'variant_importance.tsv.gz',sep='\t',index=False)
    # 原第324行：根据窗口设置枚举预期窗口数量。
    expected = sum(1 for _ in windows(v,config['window'],config['step']))
    # 原第325行：统计目前存在的窗口结果文件数量。
    done = len(list((root/'windows').glob('*/importance.csv')))
    # 原第326行：保存覆盖摘要。
    save_json(root/'coverage.json',{'expected_windows':expected,'completed_windows':done,
        # 原第327行：记录已评分位点数和总数；本版complete主要比较窗口文件数量，严格坐标/覆盖核查在后处理脚本中加强。
        'scored_variants':int(keep.sum()),'total_variants':len(v),'complete':done==expected,
        # 原第328行：明确写出仅有重要性和排名，没有校准P值、选择结论或Matrix Profile。
        'warning':'Ranks and importance only. No calibrated P values, no selection claims, no Matrix Profile region inference.'})
    # 原第329行：打印已完成/预期窗口数与已评分/全部位点数。
    print(f'Aggregated {done}/{expected} windows, {int(keep.sum())}/{len(v)} SNPs')

# 原第331行：定义命令行总入口。
def main():
    # 原第332行：创建参数解析器，帮助说明来自文件顶部docstring。
    parser = argparse.ArgumentParser(description=__doc__)
    # 原第333行：要求输入一个子命令；子命令名存到a.command。
    sub = parser.add_subparsers(dest='command',required=True)
    # 原第334行：注册prepare子命令。
    p = sub.add_parser('prepare')
    # 原第335行：prepare必须提供环境文件、VCF和输出目录。
    p.add_argument('--env',required=True);p.add_argument('--vcf',required=True);p.add_argument('--out',required=True)
    # 原第336行：可选限制过滤前读取记录数；0为全部。
    p.add_argument('--max-records',type=int,default=0)
    # 原第337行：注册structure子命令，要求准备数据目录和背景VCF。
    p = sub.add_parser('structure');p.add_argument('--data',required=True);p.add_argument('--vcf',required=True)
    # 原第338行：批量注册tune、scan、aggregate三个子命令。
    for cmd in ['tune','scan','aggregate']:
        # 原第339行：这三个都需要--data和--out。
        p = sub.add_parser(cmd);p.add_argument('--data',required=True);p.add_argument('--out',required=True)
        # 原第340行：只有调参和扫描需要训练参数，汇总不需要。
        if cmd != 'aggregate':
            # 原第341行：定义窗口宽度默认1000、训练上限默认100。
            p.add_argument('--window',type=int,default=1000);p.add_argument('--epochs',type=int,default=100)
            # 原第342行：定义默认随机种子。
            p.add_argument('--seed',type=int,default=20261023)
        # 原第343行：仅tune增加--trials，默认30。
        if cmd == 'tune':p.add_argument('--trials',type=int,default=30)
        # 原第344行：以下参数仅属于scan。
        if cmd == 'scan':
            # 原第345行：默认步长500、每窗口重复训练5次。
            p.add_argument('--step',type=int,default=500);p.add_argument('--repeats',type=int,default=5)
            # 原第346行：默认不限制窗口数量；limit-windows=2时只遍历前两个。
            p.add_argument('--limit-windows',type=int,default=0)
    # 原第347行：解析终端命令，得到参数对象a。
    a = parser.parse_args()
    # 原第348行：逐项检查部分数值参数。
    for name in ['window','epochs','trials','repeats']:
        # 原第349行：若参数存在且小于1，则报错退出；step的合法范围在windows中另外检查。
        if hasattr(a,name) and getattr(a,name)<1:parser.error(f'{name} must be positive')
    # 原第350行：根据a.command调用同名Python函数；例如tune调用tune(a)，而不是直接在命令行执行官方函数。
    globals()[a.command](a)

# 原第352行：只有直接作为脚本运行时才进入这里；被import时不自动执行main。
if __name__ == '__main__':
    # 原第353行：启动命令行解析和相应处理函数。
    main()
