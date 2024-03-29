# 背景

## Intro
- 在网络测量中，之前以理论为指导的的近似算法为减小误差需要足够的资源，于是资源配置和精准性参数有紧密联系，这导致和实践上的限制：
  - 管理员需要指定可容忍错误级别作为资源配置的输入
  - 每个配置都跟特定的参数相关联(譬如阈值)
  - 理论分析仅提供最坏情况分析，无法根据实际工作量指导资源配置;
  - 对于特定的流定义，每个配置都是固定的
  - 管理员不能轻易地量化度量结果的错误程度。
- 对于sketch算法，论文提供了一个新的测量框架以解决这种局限性
  - 它的思想是在sketch中描述资源冲突的固有统计特性，而不是追求一个完美的资源配置来缓解资源冲突

## 设计需求
- R1: 占用内存小
  - 高内存占用会加剧芯片占用和热消耗，从而增加制造成本。且即是是现代的设备的SRAM大小仍限制每流追踪
  - 服务器中的软件交换机可以利用丰富的服务器端DRAM。然而，高内存使用量不仅会耗尽共存应用程序(例如vm或容器)的内存，而且还会由于更严重的缓存争用而降低它们的性能。
- R2: 数据包处理速度快
  - 一条充分利用的10Gbps链路对应64字节数据包的速率为14.88Mpps;同样，在3GHz处理器中，每个包的时间预算只有大约200个CPU周期。过慢会导致丢包
- R3: 实时响应
  - 一些异常检测之类的任务需要快速反馈
- R4: 通用性
  - 为每种类型的流量统计设计和部署特定的解决方案是低效的，这还导致需要在不同的解决方案之间进行复杂的资源分配，以提供准确性保证。
  - 一个实用的测量解决方案应该适用于一般类型的统计。
- R5: 此外，硬件和软件之间的资源需求也有所不同。例如，ASIC交换机具有高吞吐量但内存有限;相比之下，普通主机通常有足够的内存，但CPU处理能力有限。实现软件和硬件平台的统一度量是至关重要的。

## 现存方法的问题
- L1: 很难指定“预期的误差”
  - 理论上总是证明误差小于某个特定值的可能性是多少，但是实践中如何选定这个误差让性能达到最好需要不同场景的领域知识。
- L2: 很难query不同的阈值
  - 很多需要检测关于特定阈值的统计量的问题所实现的配置都跟阈值的具体值强相关，这样导致理论分析和实现的精度都与阈值的选取相关性很大，通常在阈值过低的时候原来的方法就不好使了
- L3: 很难通过理论调整配置
  - 理论成果的应用仍然很有挑战性
    - 一些测量方法只提供渐近复杂性结果，而不是配置的封闭形式参数。
    - 其次，大多数近似测量方法执行最坏情况分析，并且不考虑实际的工作量。
    - 一些工作需要精确的分布模型输入，不能容易地适应实际的网络状态
    - SCREAM根据之前的运行时行为动态调整配置，但它需要一个集中的控制器进行复杂的协调来收集sketch统计信息和微调sketch配置。
  - 论文对比了完全使用理论进行配置和依靠经验进行配置两种情况，结果是实际的内存和CPU使用量比理论提出的要少得多。然而，调优资源高效配置是一项劳动密集型的工作。
- L4: 很难re-define flowkeys
  - 如果传到设备中的流同时又被5元组和源-目的地标识的流，那么不得不开两个配置，尽管他们使用的是相同的算法
  - 一些分层Heavy hitter检测方法可以检测多个级别的流量键，其中一个级别是另一个级别的前缀，但这些级别必须预先定义。
- L5: 很难验证正确性
  - 理论通常证明譬如“heavy hitter的算错率很小”，但不能验证“一个特定的流是否属于heavy hitter的概率很小”

# Sketchlearn

## 简介与特性
- **简介**
  - 之前的算法都根据精确度需求分配大小，而Sketchlearn分配给每个sketch固定大小的matrix of counters，并且附带hash packet或记录每一行计数器个数的字节
  - 现有的基于sketch的测量方法关注的是如何预先分配所需的最小sketch尺寸，以满足精度要求。相反，SketchLearn采用了一种全新的方法，通过统计建模来描述和过滤哈希冲突的影响。它不需要针对特定的测量任务或需求(例如，预期错误、查询阈值或流定义)微调其配置。它是自适应的，通过统计建模，通过单一配置来完成各种测量任务和需求。
- **特性**
  - 为追踪每一位bit设置多级sketch
    - 维护一个由多个sketch组成的多层次sketch，对于每个流键，每个sketch跟踪一个bit的统计量，这样，结合统计建模，sketchlearn减少了sketch的大小和需求的资源(前三个需求)，并且可以灵活地适应流键的定义(第四个局限)
  - 大小流分离
    - 多级sketch提供了一个关键属性，如果没有大流量，它的计数器值遵循高斯分布。基于此特性，SketchLearn从多级sketch中提取大流量，并将剩余计数器留给小流量，形成高斯分布。这种分离使SketchLearn能够解决各种流量统计的哈希冲突(R4)。例如，SketchLearn只将提取的大流量用于Heavy Hitter检测，但在估计基数时包括剩余计数器。
    - 这和其他测量工具不同的是，其他的测量工具设计的是为了不同的环境(譬如测量top-k)
  - 不使用参数进行模型推断
    - SketchLearn自动挖掘隐藏在多层sketch中的信息，在不依赖预期的错误值、阈值参数或调整配置(前三个之前工作的局限)的情况下区分大的和小的流。它迭代学习多层sketch内的统计分布，并利用分布来指导大流量提取。
    - 注意，无参数模型推断并不意味着SketchLearn本身是无配置的。SketchLearn仍然需要管理员配置多级sketch。此外，它不能消除查询引起的参数(Heavy-hitter或阈值)。然而，SketchLearn最小化了配置的影响，因为它的配置现在可以通过自适应建模很容易地参数化。
  - 为指定流附加错误测量
    - 在模型推断过程中，SketchLearn进一步将每个流与与给定类型的流量统计相对应的错误度量相关联。因此，我们不需要指定任何错误(L1已解决);相反，我们使用附加的错误来量化流的正确性(L5已解决)。
  - 可以基于整个网络进行推断
    - SketchLearn可以方便地总结部署在多个测量点的一个或多个多层次sketch的结果，对整个网络形成无参数的模型推断。

## 结构设计
- 网络结构
  - SketchLearn由分布式数据平面和集中控制平面组成，方式和软件定义网络测量结构类似
  - 数据平面由跨网络的多个测量点(例如，软件/硬件交换机或终端主机)组成。它在每个测点部署一个多级sketch，对进入的数据包进行处理，将数据包的统计信息记录到多级sketch中，并将多级sketch报告给控制平面进行分析。
  - 控制平面分析每个sketch并把他们分解成3个部分
    - large flow list，它指明一组大的流量，并包括每个被指明大流的估计频率和相应的误差测量
    - 残留的多层sketch，他存储过滤掉被确定的大流后剩下小流的统计量
    - bit-level的计数器值的分布，每个分布为残留的多层sketch的每一bit的分布建模
    - 控制平面有一个query runtime，它接受常规的或特别的测量查询，并为所求的特定测量任务报告来自一个或多个多层sketch(取决于查询)的相关统计信息。
  - 下图是网络结构的一个例图
  <img src=".\Sketchlearn_figure3.png" width="100%" height="100%">
- 多层sketch结构
  - 对于长为l bit的flowkey，多层sketch维持一个l+1级sketch，从第0级到第l级，每级sketch大小都为$r\times c$，level-0记录所有包的统计值，level-k($1\leq k\leq l$)记录flowkey第k bit的统计值，$V_{i,j}[k]$表示第k级sketch第i行第j列的统计值，所有sketch共享r个哈希函数$h_1,h_2,...,h_r$(注意，大多数基于sketch的方法也假定使用成对独立的哈希函数)
  - 对于key-value(f,v)，我们更新level-0 sketch的r个counter，并且当且仅当f的第k个bit为1的时更新第k个sketch
  - 下图是一个更新的例子
  <img src=".\Sketchlearn_figure4.png" width="100%" height="100%">

## Model Learning

### 动机和符号
- 我们的主要目标是建立一个统计模型，它只利用嵌入在多层sketch中的信息来减少错误，而不采取任何额外的信息。我们的模型组件考虑了两个因素。
  - 哈希冲突
    - 误差主要由哈希冲突引起，为了刻画哈希冲突的影响，我们首先为每个counter所持有的流量数目建模
    - 记号：
      - 具有相同行号和列号的counter共享相同的冲突流，我们称叫这l+1个counter一个堆栈，并用(i(行号),j(列号))表示这个堆栈
      - 用$n_{i,j}$表示哈希到栈(i,j)的流的数目，n是总的观测流数
  - bit-level的flowkey分布
    - 流键通常在其位模式中表现出非均匀的概率分布。例如，数据中心中的IP地址往往共享相同的前缀;协议字段可能等于6，因为TCP流量占主导地位。我们在bit的粒度上描述流键分布。
    - 记号
      - $f[k]$表示流f的第k位，$p[k]$表示流键的第k位为1的可能性，$p_{i,j}[k]$表示哈希到栈(i,j)的flowkey第k bit为1的概率
  - 分析方法
    - 我们的分析集中处理了上述两个因素，而不是单独处理
    - 我们用$R_{i,j}[k]=\frac{V_{i,j}[k]}{V_{i,j}[0]}$表示$V_{i,j}[k]$所记录的频率覆盖所有哈希到栈(i,j)频率的比例
    - 用$s_f$表示一个流真实的频率，于是有$V_{i,j}[k]=\sum_{h_i(f)=j}f[k]s_f$
  - 我们强调，上面定义的符号只是为了便于分析而引入的。不需要手动提前将它们的值参数化。

### Theory
- 假设
  - 我们假设哈希函数是随机的，于是$n_{i,j}\thickapprox \frac{n}{c}\ for\ 1\leq j\leq c$
  - 我们假设$h_i(f)$和$p[k]$是几乎不相关的，于是有$p_{i,j}[k]=p[k]$
- THM1：对于栈(i,j)和层级k，如果栈(i,j)没有被哈希进大流，那么$R_{i,j}[k]$服从高斯分布$N(p[k],\sigma_{j}^2[k])$，其中$\sigma_{j}^2[k]=p[k](1-p[k])/n_{i,j}$，于是，总的流服从$N(p[k],\sigma^2[k])$，其中$\sigma^2[k]=p[k](1-p[k])c/n$
  - 这是因为，如果用$X_{i,j}[k]$表示哈希到stack (i,j)的第k位的频率的随机变量，那么有$X_{i,j}[k]=\sum_{h_i(f)=j}f[k]$，而f[k]服从伯努利分布，$X_{i,j}[k]$服从二项分布，其均值为$n_{i,j}p[k]$，方差为$n_{i,k}p[k](1-p[k])$，故当$n_{i,j}$足够大的时候，$X_{i,j}[k]\thicksim N(n_{i,j}p[k],n_{i,j}p[k](1-p[k]))$(中心极限定理)，故$X_{i,j}[k]/n_{i,j}\thicksim N(p[k],p[k](1-p[k])/n_{i,j})$

### 模型推断：根据sketch得到的信息做进一步推理
- 上述定理需要哈希中没有大流，所以我们希望在多层sketch中捕获大流，为进一步刻画，我们提出三个概念
  - Large flow list F：它指明一组大的流量，并包括每个被指明大流的估计频率和相应的误差测量
  - 残留的多层sketch，他存储过滤掉被确定的大流后剩下小流的统计量
  - l个bit-level的计数器的分布$\{N(p[k],\sigma^2[k])\}$，描述$R_{i,j}[k]$的分布，被均值$p[k]$和方差$\sigma^2[k]$表征

- **算法框架**
<img src=".\Sketchlearn_Algorithm1.png" width="100%" height="100%">
  - 输入为整个多级sketch，$S=\{V_{i,j}[k]\}$
  - 首先初始化大流列表为空，控制参数$\theta = 1/2$，并且计算l个bit-level的计数器值分布(尽管有可能存在大流)(前3行)
  - 为了消除来自大流的干扰，SketchLearn基于(不准确的)分布迭代地提取大流(5-8行)
  - 在每次迭代中，它删除集合F'中提取的流，并重新计算分布(9-10行)
  - 迭代停止当且仅当新的分布很好地符合高斯分布(11-12行)
  - 如果没有捕获到大流，那么把$\theta$减半(13-14行)
  - 这样，我们需要具体地实现4个功能
    - 计算位级计数器分布(第2行和第10行)
    - 提取大流量(第7行)
    - 从S中移除提取的流量(9行)
    - 检查终止条件(11行)。

- **算法四个功能的实现**
  - 计算分布
    - 我们通过计算第k个bit统计信息的均值$p[k]$和方差$\sigma^2[k]$计算分布$N(p[k],\sigma^2[k])$
    - 注意到，对于每个i,j，$R_{i,j}[k]$都可以视为一次采样，所以我们用所有的$R_{i,j}[k]$计算关于k的均值和方差
    - 理论证明，对于高斯分布，这种估值是无偏的，方差最小
    - 如果大流都被剔除，估值将收敛于真值
  - 提取大流(如果这一部分看不懂的话，论文里figure5 举了一个例子，从第582页的左下角开始)
    - 它是模型推理的核心子程序。它在**每个堆栈**的基础上工作
    - 输入是一个阈值$\theta$，一个想要在其中找大流的栈(i,j)，多级sketch S，估计的分布
    - 直觉是一个大流显著影响S中的计数器值。当堆栈(i,j)中的其他流具有有限的大小时，一个大流支配堆栈，通常使$R_{i,j}[k]$要么非常大(如果它的第k位是1)，要么非常小(如果它的第k位是0)。即使大流量不是明显的主导，它至少应该改变计数器分布，使$R_{i,j}[k]$与它的期望p[k]有很大的偏差。因此，通过检验$R_{i,j}[k]$及其与p[k]的差值，我们可以确定大流量是否存在。
    - 整个子算法一共有5个部分
      - 第一部分(进入这个算法算法说明我们假定了假定栈(i,j)存在大流)，寻找大流在第k位bit为1的概率$\hat p[k]$
        - 如果$R_{i,j}[k]<\theta$，说明大流没有哈希进来，置$\hat p[k]=0$
        - 如果$1-R_{i,j}[k]<\theta$，说明大流大概率被哈希进来了，置$\hat p[k]=1$
        - else, 我们假设存在一个大流第k位为1(这个流的大小至少为$V_{i,j}[0]\theta$)，我们可以通过模拟驱逐掉这个大流之后剩余的$R_{i,j}[k]$和$p[k]$之间的差值使用贝叶斯理论计算出有大流哈希到此处相应的概率$\hat p[k]$
      - 第二部分，找到所有可能为大流的候选流键
        - 我们逐位确定相应大流的bit
        - 如果$\hat p[k]$趋近1(论文中使用>0.99判定)，则断言存在于这个栈中的大流这一位为1，趋近0(<0.01)则断言0
        - else, 我们把这一位置为*，以后筛选到底哪个流是大流的时候枚举这一位的两种可能
        - 得到所有候选流之后，我们check其中的流是否会哈希到(i,j)(应该是要剔除掉不会哈希到此处的流)
      - 第三部分，估计大流的频率
        - 我们根据“剔除掉这个流之后剩下的R值应该趋近于p”来估算最可能的大流频率
        - 我们希望在把大流扔出去之后，R变为p，所以得到以下估计
        - $$s_f=\begin{cases}
            \frac{R_{i,j}[k]-p[k]}{1-p[k]}V_{i,j}[0] & if\ f[k]=1\\
            (1-\frac{R_{i,j}[k]}{p[k]})V_{i,j}[0] & if\ f[k]=0
        \end{cases}$$
        - 最后我们取估值为所有l次估值后的中值作为最终估计
      - 第四部分，把得到的flowkey和一个bit-level的可能性向量联系起来
        - 这个向量的第k位是$\hat p[k]$(如果选择了这一位为1)或$1-\hat p[k]$(如果没选择)
        - 这个向量中的值越大，意味者这个流是真的大流的概率越大
        - 这个向量在之后可以被用来检测错误
      - 第五部分，验证候选的大流
        - 因为有可能匹配*的时候匹配到不是大流的大流，最后我们再筛一遍拿到的所有可能是大流的大流
        - 具体的做法是把对应的流去其他不同行的栈中哈希一遍，选择当前估值和其他栈中哈希得到的最小的估值作为最终估值，当且仅当最终估值小于$\theta V_{i,j}[0]$的时候我们才声称捕获了这个大流
    - 下面是一个例图
    <img src=".\Sketchlearn_figure5.png" width="100%" height="100%">
  - 从S中删去大流
    - 把哈希到的栈中f[k]=1的计数器中的值减去f的估计值
    - 由于检测出假的大流的可能被上一步的re-hash减少了，所以这一步进行错误的删除的可能性很小，对最终结果的影响也很小
  - 检查终条件
    - 我们要检测现在是否存在大流，即检测现在的分布是否已经符合高斯分布了
    - 按照高斯分布，当$68.26\%,95.44\%,99.73\%$的$R_{i,j}[k]$的观测值的差距在1，2，3阶标准差之内的时候，我们终止算法

- **算法分析**
  - 请注意，当只有少量流时，每个堆栈中很少有哈希冲突。在这种情况下，每个流很可能支配它自己的哈希堆栈，可以很容易地提取。因此，我们的分析只集中在有大量小流的情况下。
  - 我们的算法精确性和捕获到的大流的最小值相关，这个值越小，意味者可以直接捕获到的大流更多，并且可以让筛剩下的值越小，从而更加适应正态分布
  - 可以保证的捕获到的频率
    - THM2：对于有c列的多级sketch，有多于总流量$\frac{1}{c}$的频率的流必然会被捕获
      - 这个证明写了半页orz 具体的证明在[这里](https://github.com/huangqundl/SketchLearn/blob/master/TechReport.pdf)
    - 一个简单的版本是，对任意栈，超过总流量一半的流必然会被检测出来，这是因为一开始把$\theta$设置成了1/2，某个流超过一半意味着它bit值为1的位必然导致$1-R_{i,j}[k]<\theta$，为0的位必然导致$R_{i,j}[k]>\theta$，从而它必然进入候选大流队列，而且不会被筛出去
  - 推荐配置
    - 管理员需要为多级sketch配置行数和列数。对于行数，我们的评估表明，一行就足够了。对于列数，可以根据内存预算配置，也可以指定一个保证频率，并使用THM2将其转换为列数。
  - 内参数
    - 模型推理引入了两个内部参数。注意，管理员不需要关心它们的实际值;事实上，根据我们的评估，他们当前的设置对各种类型的流量统计都很有效
    - 第一个参数是$\theta$，在开始被设为1/2，初始值的动机是THM2，它保证了有多大的流量被提取，并使分析剩下的sketch的时候可以快速收敛。为控制捕获的过程，对后续的$\theta$值进行迭代减半。我们的经验是，我们的大流量提取程序适用于任何$\theta$递减序列。
    - 第二个参数是当$\hat p[k]$多大时我们选取这个bit的阈值(上面我们提到的0.99)，其原理是我们要给概率小于这个值的k位赋值*。较大的值会导致更多的∗，并最终减少错误的位分配的错误。
  - Discussion
    - 我们的模型解决了L1到L3的限制。
    - L1: 首先，定理1和定理2共同保证了模型推理的正确性，因此我们不需要配置任何期望的误差概率作为输入。虽然SketchLearn仍然使用两个用户指定的参数，但这两个参数都很容易配置，因为管理员可以很容易地知道他们设备中的最小流量大小或内存预算。
    - L2: 其次，对于小于总频率1/c的流量，SketchLearn也努力提取它们，尽管它们在理论上没有保证。我们的经验是，即使需要提取保证频率的50%作为需要的阈值(即1/2c)，仍有超过99%的流量被提取出来。因此，SketchLearn对非常小的阈值是robust的。
    - L3: 最后，我们的模型推理是迭代的，直到结果符合两个定理。因此，管理员只需遵循这些定理，而不需要进一步调优。

## Query Runtime
- query 每个流的频率
  - Sketchlearn首先查找大流队列，如果找到了，返回估计值，否则查找剩余的sketch，返回可靠性向量，如同提取大流的第三和第四部分
- query Heavy hitters
  - 因为我们能将足够小下限大小的大流都捕获到，因此，只需要查找大流队列有没有超过给定阈值的即可(要强调的是这个阈值跟主体算法并没有关系)
- query Heavy changers
  - 我们直接对于两个epoch都查一遍heavy hitters，然后声称其中变化超过给定阈值的即为heavy changers
- 估计流的基数(求n)
  - 我们已经得到了每个过滤大流后的sketch的分布，所以直接代入THM1中的式子，对每个k，求出对应的n的估值$n=\frac{p[k](1-p[k])c}{\sigma^2[k]}$，然后取其中间值，然后在此基础上加上所有被放在大流队列中的所有流的估值大小即可
- 频率的分布和熵
  - 频率分布和熵。残差(即剔除大流之后)的sketch压缩了计数器中流量频率分布的信息。SketchLearn使用MRAC技术恢复分布，该技术使用期望最大化来匹配观察到的计数器值。有了频率分布，许多其他统计量，如熵，也可以计算。管理员还可以使用其他评估方法。

## Extended Queries
- SketchLearn提出了两个扩展，在测量结果上附加误差测量和允许使用更多流键的定义
  - 在测量结果上附加误差测量
    - 统计量的估计数据来自残差sketch和大流队列，我们处理了他们的误差估计
    - 对于残差sketch的误差，一方面是来自对高斯分布的适应性，而因为我们保证了筛选出了大流，因此适应性是良好的，另一方面来自对均值p和方差$\sigma^2$的估计，而理论上这种采样是无偏并且方差最小的，因此，这方面的误差可以丢弃，忽略不计
    - 对于大流队列的估计误差，是由于假阳性的存在。Sketchlearn将bit-level的误差向量提供给管理者，管理员可以使用误差向量决定是否排除某些流，譬如可以排除具有一半以上概率$<90\%$的分量的流，评估表明这个简单的过滤器有效地消除了几乎所有的误报。管理员还可以根据需要应用其他过滤器。我们声称这样一个错误过滤器比在度量之前指定资源配置中的错误要容易得多(参见Limitation L1)，因为管理员可以在收集度量结果之后调整错误过滤器
  - 允许使用更多流键的定义
    - SketchLearn可以使用任何原始流键定义的流键(例如，5元组流中104位的任何子集)对其中任何比特组合进行流量统计。对于大流列表，SketchLearn可以基于流的定义简单地加起来所有值;对于后者,SketchLearn只依据感兴趣的bit位做检测
    - 一个问题是也许大流会因为不同表示而被分解成若干不同的部分，但是有的部分没有被哈希到相同的stack中，这就会导致有些流量没有被估计出来。我们说明通常情况下大流会在每个别名中都具有很大的流量，所以这个问题被避免了

## 互联网级别的合作
- SketchLearn允许管理员访问部分或所有测量点，以计算网络范围内的测量统计信息。特别是，它方便地支持全网部署和全网集成。
- 全网部署
  - SketchLearn可以部署在任何硬件/软件交换机或服务器主机上。通过跨测量点改变哈希函数，我们实际上为多级sketch提供了更多的行。因此，在实践中，在每个测量点为每个多级草图分配一行就足够了。此外，我们不限制部署决策，使SketchLearn正交于以前的网络测量部署方法。
- 全网集成
  - 因为一个流通常会遍历多个测量点，所以我们可以通过利用来自其他测量点的结果来改进我们对特定测量点的模型学习。
  - 具体地说，我们检查每个提取的流沿其遍历路径的比特级概率。如果一个流在路径中的大多数测量点上都不能通过我们的错误过滤器，我们就丢弃它。另一方面，如果一个流只在一个测量点丢失，而它在沿着路径的其他测量点的比特级概率很高，我们查询它在miss掉的测量点的bit-level的概率向量，如果这个向量和其他测量点的一致，那我们就接收这个大流
  - 如何聚合来自多个测量点的结果取决于测量任务和网络拓扑，我们将决策留给管理员。

## 总结
- 通过解耦资源配置和精度参数之间的绑定，SketchLearn为近似测量提供了一个新的视角。它的想法是利用自动统计推断来提取流量统计信息。实验表明，SketchLearn是资源高效、准确和通用的，而且在实践上对用户的要求非常少。