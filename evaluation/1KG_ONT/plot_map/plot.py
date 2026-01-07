import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
import numpy as np

# Original data structure
all_populations = {
    "Esan in Nigeria": (0, 0, 48, "ESN"),
    "African Caribbean in Barbados": (0, 0, 46, "ACB"),
    "Colombian in Medellin, Colombia": (0, 0, 44, "CLM"),
    "Peruvian in Lima, Peru": (0, 0, 44, "PEL"),
    "South Africa": (0, 0, 43, "BEB"),
    "Kinh in Ho Chi Minh City, Vietnam": (0, 0, 43, "KHV"),
    "Puerto Rican in Puerto Rico": (0, 0, 42, "PUR"),
    "Southern Han Chinese, China": (0, 0, 41, "CHS"),
    "British in England and Scotland": (0, 0, 41, "GBR"),
    "Iberian populations in Spain": (0, 0, 40, "IBS"),
    "Sri Lankan Tamil in the UK": (0, 0, 40, "STU"),
    "Punjabi in Lahore,Pakistan": (0, 0, 40, "PJL"),
    "Indian Telugu in the UK": (0, 0, 39, "ITU"),
    "Chinese Dai in Xishuangbanna, China": (0, 0, 39, "CDX"),
    "Mende in Sierra Leone": (0, 0, 39, "MSL"),
    "Gujarati Indian in Houston,TX": (0, 0, 37, "GIH"),
    "Finnish in Finland": (0, 0, 36, "FIN"),
    "Yoruba in Ibadan, Nigeria": (0, 0, 36, "YRI"),
    "African Ancestry in Southwest US": (0, 0, 36, "ASW"),
    "Luhya in Webuye, Kenya": (0, 0, 35, "LWK"),
    "Han Chinese in Bejing, China": (0, 0, 35, "CHB"),
    "Mexican": (0, 0, 34, "MXL"),
    "Toscani in Italy": (0, 0, 34, "TSI"),
    "Gambian in Western Division, The Gambia": (0, 0, 33, "GWD"),
    "Japanese in Tokyo, Japan": (0, 0, 32, "JPT"),
    "Utah residents with Northern and Western European ancestry": (0, 0, 29, "CEU")
}

# Convert to DataFrame
df_all = pd.DataFrame(all_populations).T
df_all.columns = ['Longitude', 'Latitude', 'Population', 'Code']

# Ensure 'Population' is numeric
df_all['Population'] = pd.to_numeric(df_all['Population'], errors='coerce')

# Create a mapping of population names to countries
country_mapping = {
    "Nigeria": ["Esan in Nigeria", "Yoruba in Ibadan, Nigeria"],
    "Barbados": ["African Caribbean in Barbados"],
    "Colombia": ["Colombian in Medellin, Colombia"],
    "Peru": ["Peruvian in Lima, Peru"],
    "South Africa": ["South Africa"],
    "Vietnam": ["Kinh in Ho Chi Minh City, Vietnam"],
    "Puerto Rico": ["Puerto Rican in Puerto Rico"],
    "China": ["Southern Han Chinese, China", "Chinese Dai in Xishuangbanna, China", "Han Chinese in Bejing, China"],
    "United Kingdom": ["British in England and Scotland", "Sri Lankan Tamil in the UK", "Indian Telugu in the UK"],
    "Spain": ["Iberian populations in Spain"],
    "Pakistan": ["Punjabi in Lahore,Pakistan"],
    "Sierra Leone": ["Mende in Sierra Leone"],
    "United States": ["Gujarati Indian in Houston,TX", "African Ancestry in Southwest US", "Utah residents with Northern and Western European ancestry"],
    "Finland": ["Finnish in Finland"],
    "Kenya": ["Luhya in Webuye, Kenya"],
    "Mexico": ["Mexican"],
    "Italy": ["Toscani in Italy"],
    "Gambia": ["Gambian in Western Division, The Gambia"],
    "Japan": ["Japanese in Tokyo, Japan"]
}

# Create a new DataFrame for countries
df_countries = pd.DataFrame(columns=['Population', 'Code'])

for country, populations in country_mapping.items():
    total_population = df_all.loc[df_all.index.isin(populations), 'Population'].sum()
    country_code = df_all.loc[df_all.index.isin(populations), 'Code'].iloc[0]  # Use the first code found
    df_countries.loc[country] = [total_population, country_code]

# Superpopulation mappings
population_superpopulation_dict = {
    'CDX': 'EAS',  # Chinese Dai in Xishuangbanna, China
    'CHB': 'EAS',  # Han Chinese in Bejing, China
    'JPT': 'EAS',  # Japanese in Tokyo, Japan
    'KHV': 'EAS',  # Kinh in Ho Chi Minh City, Vietnam
    'CHS': 'EAS',  # Southern Han Chinese, China
    'BEB': 'SAS',  # Bengali in Bangladesh
    'GIH': 'SAS',  # Gujarati Indian in Houston,TX
    'ITU': 'SAS',  # Indian Telugu in the UK
    'PJL': 'SAS',  # Punjabi in Lahore,Pakistan
    'STU': 'SAS',  # Sri Lankan Tamil in the UK
    'ASW': 'AFR',  # African Ancestry in Southwest US
    'ACB': 'AFR',  # African Caribbean in Barbados
    'ESN': 'AFR',  # Esan in Nigeria
    'GWD': 'AFR',  # Gambian in Western Division, The Gambia
    'LWK': 'AFR',  # Luhya in Webuye, Kenya
    'MSL': 'AFR',  # Mende in Sierra Leone
    'YRI': 'AFR',  # Yoruba in Ibadan, Nigeria
    'GBR': 'EUR',  # British in England and Scotland
    'FIN': 'EUR',  # Finnish in Finland
    'IBS': 'EUR',  # Iberian populations in Spain
    'TSI': 'EUR',  # Toscani in Italy
    'CEU': 'EUR',  # Utah residents with Northern and Western European ancestry
    'CLM': 'AMR',  # Colombian in Medellin, Colombia
    'MXL': 'AMR',  # Mexican Ancestry in Los Angeles, California
    'PEL': 'AMR',  # Peruvian in Lima, Peru
    'PUR': 'AMR'   # Puerto Rican in Puerto Rico
}

# Map superpopulation to the DataFrame
df_all['Superpopulation'] = df_all['Code'].map(population_superpopulation_dict)

# Create custom colormap
colors = ["#C6EDCF", "#22C9BF", "#49C7DA", "#006BB2", "#002786"]
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)

# Load shapefile
shapefile_path = "ne_110m_admin_0_countries.shp"
world = gpd.read_file(shapefile_path)

# Function to assign population values to countries
def assign_population(x):
    country_name = x['NAME']
    if country_name in df_countries.index:
        return df_countries.loc[country_name, 'Population']
    return None

# Assign population values to countries
world['Population'] = world.apply(assign_population, axis=1)

# Filter out countries with no data
world = world.dropna(subset=['Population'])

# Create figure with adjusted size
fig = plt.figure(figsize=(27, 12))

# Define the width ratio (4:1) with increased gap
total_width = 0.82
map_width = total_width * 4/5
bar_width = total_width * 1/5
gap = 0.08

# Create map subplot
ax_map = fig.add_axes([0.05, 0.1, map_width, 0.8], projection=ccrs.PlateCarree())

# Draw world map
world.boundary.plot(ax=ax_map, linewidth=1)

# Fill country areas
world.plot(column='Population', cmap=custom_cmap, linewidth=0.8, ax=ax_map, edgecolor='0.8')

# Add coastlines and borders
ax_map.add_feature(cfeature.COASTLINE, linewidth=0.3)
ax_map.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.2)

# Set map extent
ax_map.set_extent([-180, 180, -60, 90], crs=ccrs.PlateCarree())

# Add x and y ticks
ax_map.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
ax_map.set_yticks([-60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
ax_map.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, p: f'{v:.0f}°'))
ax_map.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, p: f'{v:.0f}°'))

# ============ 调整地图刻度字体大小 ============
ax_map.tick_params(axis='both', labelsize=25)  # 可以根据需要调整数值

# Set labels for world map
ax_map.set_xlabel('Longitude', fontsize=25, fontweight='bold')
ax_map.set_ylabel('Latitude', fontsize=25, fontweight='bold')
ax_map.set_title('Global Population Distribution', fontsize=25, fontweight='bold')

# Add manual scale bar
scale_bar_length = 2000  # km
scale_bar_x = -170
scale_bar_y = -50
ax_map.plot([scale_bar_x, scale_bar_x + scale_bar_length/111], [scale_bar_y, scale_bar_y], 
            color='black', linewidth=3, transform=ccrs.PlateCarree())
ax_map.text(scale_bar_x + scale_bar_length/222, scale_bar_y - 3, f'{scale_bar_length} km', 
            horizontalalignment='center', verticalalignment='top', 
            transform=ccrs.PlateCarree(), fontsize=20, fontweight='bold')

# Add colorbar
sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=Normalize(vmin=world['Population'].min(), vmax=world['Population'].max()))
sm._A = []
cbar = fig.colorbar(sm, ax=ax_map, orientation='horizontal', pad=0.08, aspect=30)
cbar.set_label('Population Count', fontsize=22, fontweight='bold')
# ============ 调整 colorbar 刻度字体大小 ============
cbar.ax.tick_params(labelsize=22)  # 可以根据需要调整数值

# ======================= Modified Bar Plot Section =========================

# Create bar plot subplot with increased left position (moved right)
ax_bar = fig.add_axes([0.05 + map_width + gap, 0.1, bar_width, 0.8])

# Sort data by superpopulation and population
sorted_df = df_all.sort_values(['Superpopulation', 'Population'], ascending=[True, False])

# Reset index to ensure proper ordering
sorted_df = sorted_df.reset_index()

# Define superpopulation categories and their corresponding colors for background
superpop_categories = sorted_df['Superpopulation'].unique()
superpop_colors = {
    'EAS': '#8f96bd',
    'SAS': '#d6d69a',
    'AFR': '#d9e6eb',
    'EUR': '#2A347A',
    'AMR': '#9fc3d5'
}

# Calculate maximum population for setting label positions
max_population = sorted_df['Population'].max()

# Determine the y positions for each superpopulation group
current_y = 0
group_boundaries = {}
for category in superpop_categories:
    group_size = sorted_df[sorted_df['Superpopulation'] == category].shape[0]
    group_boundaries[category] = (current_y, current_y + group_size)
    current_y += group_size

# Add background rectangles for each superpopulation group
for category, (start, end) in group_boundaries.items():
    ax_bar.axhspan(start - 0.5, end - 0.5, facecolor=superpop_colors.get(category, '#ffffff'), alpha=0.6)
    ax_bar.text(-max_population * 0.55, (start + end) / 2 - 0.5, category, 
                va='center', ha='right', fontsize=22, fontweight='bold')

# Plot bars
bars = ax_bar.barh(
    y=range(len(sorted_df)),
    width=sorted_df['Population'],
    color=custom_cmap(Normalize(vmin=df_all['Population'].min(), 
                                vmax=df_all['Population'].max())(sorted_df['Population'])),
    height=0.8
)

# Set y-ticks to be in the middle of each bar
ax_bar.set_yticks(range(len(sorted_df)))
ax_bar.set_yticklabels(sorted_df['Code'], fontsize=24)

# Set bar plot labels
ax_bar.set_xlabel('Population Count', fontsize=26, fontweight='bold')
ax_bar.set_ylabel('')

# ============ 调整柱状图刻度字体大小 ============
ax_bar.tick_params(axis='x', labelsize=20)  # x轴刻度字体大小
ax_bar.tick_params(axis='y', labelsize=24)  # y轴刻度字体大小（已在上面的set_yticklabels设置）

# Remove top and right spines
ax_bar.spines['top'].set_visible(False)
ax_bar.spines['right'].set_visible(False)

# Adjust x-axis to include more space for superpopulation labels
ax_bar.set_xlim(-max_population * 0.25, sorted_df['Population'].max() * 1.1)

# Add value labels to the end of each bar
for i, bar in enumerate(bars):
    width = bar.get_width()
    ax_bar.text(width, bar.get_y() + bar.get_height()/2, f'{width:.0f}', 
                ha='left', va='center', fontweight='bold', fontsize=22)

# Hide the left spine as labels are on the far left
ax_bar.spines['left'].set_visible(False)

# Improve layout by adjusting the subplot parameters
plt.tight_layout(rect=[0, 0, 1, 1])

# Save figure as svg
plt.savefig('population_map_and_chart.pdf', format='pdf', dpi=600)

# Show plot
plt.show()