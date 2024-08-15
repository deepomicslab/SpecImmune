import numpy as np
import pandas as pd
from plotnine import ggplot, aes, geom_line, labs, theme_minimal, theme, scale_y_continuous, element_text

# 模拟数据
np.random.seed(42)
test_counts = np.arange(1, 101)
accuracy_software_1 = np.random.normal(loc=0.8, scale=0.05, size=len(test_counts))
accuracy_software_2 = np.random.normal(loc=0.75, scale=0.05, size=len(test_counts))
accuracy_software_3 = np.random.normal(loc=0.7, scale=0.05, size=len(test_counts))

# 创建数据框
data = pd.DataFrame({
    'Test_Counts': np.tile(test_counts, 3),
    'Accuracy': np.concatenate([accuracy_software_1, accuracy_software_2, accuracy_software_3]),
    'Software': np.repeat(['Software 1', 'Software 2', 'Software 3'], len(test_counts))
})

# 绘制折线图
p = (ggplot(data, aes(x='Test_Counts', y='Accuracy', color='Software'))
     + geom_line(size=1.2)
     + labs(title='Accuracy of Three Software over Multiple Tests',
            x='Test Counts',
            y='Accuracy')
     + theme_minimal()
     + theme(
         plot_title=element_text(size=16, face='bold', margin={'b': 10}),
         axis_title_x=element_text(size=14, margin={'t': 10}),
         axis_title_y=element_text(size=14, margin={'r': 10}),
         legend_title=element_text(size=12),
         legend_position='right'
     )
     + scale_y_continuous(limits=[0, 1])
    )

print(p)