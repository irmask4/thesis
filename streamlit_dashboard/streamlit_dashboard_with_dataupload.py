import streamlit as st
import pandas as pd
import networkx as nx
from pyvis.network import Network
import streamlit.components.v1 as components
from collections import defaultdict
import plotly.graph_objects as go
import seaborn as sns
from pandas.api.types import is_numeric_dtype, is_string_dtype
import plotly.express as px
import re

# --- Page setup ---
st.set_page_config(layout="wide")
st.title("GRN & Enrichment Dashboard")

# --- File upload ---
st.sidebar.header("Upload Data Files")
st.sidebar.markdown("**Important:** Please upload each file in the correct box below - the GRN file must go into the GRN uploader and the enrichment file into the Enrichment uploader. Swapping them will cause errors.")


grn_file = st.sidebar.file_uploader("Upload preprocessed GRN file:", type="csv")
enrich_file = st.sidebar.file_uploader("Upload preprocessed enrichment file:", type="csv")


# --- Load data ---
@st.cache_data
def load_data(grn_file, enrich_file):
    grn = pd.read_csv(grn_file) if grn_file else pd.DataFrame()
    enrich = pd.read_csv(enrich_file) if enrich_file else pd.DataFrame()
    return grn, enrich

grn_df, enrich_df = load_data(grn_file, enrich_file)

def data_check_grn(grn_df):
    if not grn_df.empty:
        if not (is_numeric_dtype(grn_df['Module']) or is_string_dtype(grn_df['Module'])):
            st.error("GRN 'Module' column must be either numeric or string.")
            st.stop()
        if not is_string_dtype(grn_df['Regulator']):
            st.error("GRN 'Regulator' column must be a string.")
            st.stop()
        if not is_string_dtype(grn_df['TargetGene']):
            st.error("GRN 'TargetGene' column must be a string.")
            st.stop()
        if not is_numeric_dtype(grn_df['Weight']):
            st.error("GRN 'Weight' column must be numeric.")
            st.stop()

        # Missing value checks
        for col in ['Module', 'Regulator', 'TargetGene', 'Weight']:
            if grn_df[col].isnull().any():
                st.error(f"GRN '{col}' column contains missing values. ")
                st.stop()
        for col in ['RegulatorBiotype', 'TargetBiotype']:
            if col in grn_df.columns and grn_df[col].isnull().any():
                st.error(f"GRN '{col}' column contains missing values. ")
                st.stop()

def data_check_enrich(enrich_df):
    if not enrich_df.empty:
        if not (is_numeric_dtype(enrich_df['Module']) or is_string_dtype(enrich_df['Module'])):
            st.error("Enrichment 'Module' column must be either numeric or string.")
            st.stop()
        if not is_string_dtype(enrich_df['Gene']):
            st.error("Enrichment 'Gene' column must be a string.")
            st.stop()
        if not is_string_dtype(enrich_df['Path_ID']):
            st.error("Enrichment 'Path_ID' column must be a string.")
            st.stop()
        if not is_string_dtype(enrich_df['Path_Description']):
            st.error("Enrichment 'Path_Description' column must be a string.")
            st.stop()
        if not is_string_dtype(enrich_df['Source']):
            st.error("Enrichment 'Source' column must be a string.")
            st.stop()
        # Missing value checks
        for col in ['Module', 'Gene', 'Path_ID', 'Path_Description', 'Source']:
            if enrich_df[col].isnull().any():
                st.error(f"Enrichment '{col}' column contains missing values. ")
                st.stop()    

# --- Determine views ---
views = []
if not grn_df.empty:
    views.append("GRN View")
    required_grn_cols = ['Module', 'Regulator', 'TargetGene', 'Weight']
    if not all(col in grn_df.columns for col in required_grn_cols):
        st.error(f"GRN file is missing required columns: {required_grn_cols}")
        st.stop()
if not enrich_df.empty:
    views.append("Enrichment View")
    required_enrich_cols = ['Module', 'Gene', 'Path_ID', 'Path_Description', 'Source']
    if not all(col in enrich_df.columns for col in required_enrich_cols):
        st.error(f"Enrichment file is missing required columns: {required_enrich_cols}")
        st.stop()
if not views:
    st.warning("Please upload at least one dataset to begin.")
    st.stop()

# --- Data type checks ---
data_check_grn(grn_df)
data_check_enrich(enrich_df)


# --- Sidebar navigation ---
st.sidebar.title("Navigation")
selected_view = st.sidebar.radio("Choose View", views)

# --- Module selection ---
def natural_sort(values):
    def sort_key(val):
        match = re.search(r"\d+", str(val))
        return int(match.group()) if match else float('inf')
    return sorted(values, key=sort_key)

# --- Module selection based on selected view ---
def get_available_modules(view, grn_df, enrich_df):
    if view == "GRN View" and not grn_df.empty:
        return natural_sort(grn_df['Module'].unique())
    elif view == "Enrichment View" and not enrich_df.empty:
        return natural_sort(enrich_df['Module'].unique())
    else:
        return []

module = st.sidebar.selectbox("Select Module", get_available_modules(selected_view, grn_df, enrich_df))
# --- GRN View ---
if selected_view == "GRN View":
    st.header("Gene Regulatory Network Viewer")

    threshold = st.sidebar.slider("Filter edges with |Weight| >=", 0.0, 1.0, 0.0, 0.05)
    if not enrich_df.empty:
        show_enriched_only = st.sidebar.checkbox("Show only enriched genes")
    else:
        show_enriched_only = False 

    # --- Enrichment logic ---
    enriched_map = defaultdict(list)
    enriched_genes = set()
    if not enrich_df.empty and not enrich_df[enrich_df['Module'] == module].empty:
        source_type = st.sidebar.selectbox("Enrichment Source", enrich_df['Source'].unique())
        enrich_sub = enrich_df[(enrich_df['Module'] == module) & (enrich_df['Source'] == source_type)]
        for _, row in enrich_sub.iterrows():
            enriched_map[row['Gene']].append((row['Path_ID'], row['Path_Description']))
        enriched_genes = set(enriched_map.keys())

    # --- Focus input ---
    if "focus_node" not in st.session_state:
        st.session_state.focus_node = ""

    if not show_enriched_only:
        focus_input = st.sidebar.text_input("Enter gene to focus", st.session_state.focus_node)
        if st.sidebar.button("Apply focus"):
            st.session_state.focus_node = focus_input.strip()

        focus_node = st.session_state.focus_node.strip()
    else:
        focus_node = ""

    # --- Filter GRN ---
    grn_filtered = grn_df[(grn_df['Module'] == module) & (grn_df['Weight'].abs() >= threshold)]

    # --- Build graph ---
    G = nx.DiGraph()
    for _, row in grn_filtered.iterrows():
        G.add_edge(row['Regulator'], row['TargetGene'], weight=row['Weight'])

    # Check for duplicate interactions in the selected module
    duplicates = grn_filtered.duplicated(subset=['Regulator', 'TargetGene'], keep=False)
    if duplicates.any():
        st.warning(
            "Warning: Duplicate interactions detected in this module. "
            "DiGraph will only keep one edge per interaction pair in the same direction. Check your data for duplicates."
        )

    # --- Biotype logic ---
    has_biotype = 'RegulatorBiotype' in grn_df.columns and 'TargetBiotype' in grn_df.columns
    biotype_dict = {}
    if has_biotype:
        for _, row in grn_filtered.iterrows():
            if pd.notna(row['Regulator']) and pd.notna(row['RegulatorBiotype']):
                biotype_dict[row['Regulator']] = row['RegulatorBiotype']
            if pd.notna(row['TargetGene']) and pd.notna(row['TargetBiotype']):
                biotype_dict[row['TargetGene']] = row['TargetBiotype']

    # --- Focus logic ---
    successors, predecessors, reciprocals = set(), set(), set()
    connected_nodes = set()
    if focus_node and focus_node in G:
        successors = set(G.successors(focus_node))
        predecessors = set(G.predecessors(focus_node))
        connected_nodes = successors | predecessors | {focus_node}
        for node in successors:
            if G.has_edge(node, focus_node):
                reciprocals.add(node)
        for node in predecessors:
            if G.has_edge(focus_node, node):
                reciprocals.add(node)
    # --- Handle missing focus ---
    if focus_node and focus_node not in G:
        st.warning(f"Gene '{focus_node}' is not present in selected module '{module}'.")
        st.stop()

    # --- Node selection ---
    if focus_node and focus_node in G:
        subgraph_nodes = connected_nodes
    elif show_enriched_only:
        subgraph_nodes = {node for node in enriched_genes if node in G}
    else:
        subgraph_nodes = G.nodes()

    # --- Pyvis network ---
    net = Network(height="800px", width="100%", directed=True)
    
        # Show physics buttons only if user selects checkbox
    show_physics_buttons = st.sidebar.checkbox("Show physics options", value=False)
    if show_physics_buttons:
        net.show_buttons(filter_=['physics'])


    degree_centrality = nx.degree_centrality(G)
    for node in subgraph_nodes:
        centrality = degree_centrality[node]
        size = 5 + 20 * centrality  # scale according to degree centrality
        in_deg = G.in_degree(node)
        out_deg = G.out_degree(node)

        tooltip = f"Gene: {node}\nIn-degree: {in_deg}\nOut-degree: {out_deg}"
        if node in enriched_map:
            enrichments = "\n".join([f"{desc} ({pid})" for pid, desc in enriched_map[node]])
            tooltip += f"\nEnriched in:\n{enrichments}"

        color = (
            "blue" if node == focus_node else
            "purple" if node in reciprocals else
            "orange" if node in successors else
            "#50C878" if node in predecessors else #emerald green
            "#D3D3D3"
        )

        shape = "dot"
        if has_biotype:
            if (
            grn_df['RegulatorBiotype'].eq('protein_coding').any() or
            grn_df['TargetBiotype'].eq('protein_coding').any()):
                biotype = biotype_dict.get(node, "")
                if biotype != "protein_coding":
                    shape = "triangle"

        net.add_node(node, label=node, color=color, size=size, shape=shape, title=tooltip)
        
    
    # --- Add edges ---
    for source, target, data in G.edges(data=True):
        if source in subgraph_nodes and target in subgraph_nodes:
            weight = data['weight']
            color = "green" if weight >= 0 else "red" 
            net.add_edge(source, target, color=color, width=3 * abs(weight))

    # --- Generate HTML string and render directly ---
    html_string = net.generate_html()
    components.html(html_string, height=1500, scrolling=True)
    

# --- Enrichment View ---
elif selected_view == "Enrichment View":
    st.header("Pathway Enrichment View")
    st.write(f"Module: **{module}**")

    def format_genes_multiline(gene_list, wrap=12):
        sorted_genes = sorted(set(gene_list))
        wrapped = []
        for i in range(0, len(sorted_genes), wrap):
            line = sorted_genes[i:i + wrap]
            wrapped.append(", ".join(line))
        return ",<br>".join(wrapped)

    for source in enrich_df['Source'].unique():
        st.subheader(f"{source} Enrichment")
        enrich_plot_df = enrich_df[(enrich_df['Module'] == module) & (enrich_df['Source'] == source)]

        if enrich_plot_df.empty:
            st.warning(f"No enrichment data found for {source}.")
            continue

        grouped = enrich_plot_df.groupby(['Path_ID', 'Path_Description'])
        count_df = grouped['Gene'].nunique().reset_index(name="Enriched_Genes")
        gene_lists = grouped['Gene'].apply(lambda genes: format_genes_multiline(genes)).reset_index(name="Genes")
        merged_df = pd.merge(count_df, gene_lists, on=['Path_ID', 'Path_Description'])
        merged_df.sort_values('Enriched_Genes', ascending=True, inplace=True)

        colors = sns.color_palette("viridis", n_colors=len(merged_df)).as_hex()
        fig = go.Figure()
        fig.add_trace(go.Bar(
            x=merged_df['Enriched_Genes'],
            y=merged_df['Path_Description'],
            orientation='h',
            marker_color=colors,
            customdata=merged_df['Genes'],
            hovertemplate="<b>%{y}</b><br>Genes: %{customdata}<br>Count: %{x}<extra></extra>"
        ))
        fig.update_layout(
            title=f"{source} Pathways - {module}",
            height=25 * len(merged_df) + 100,
            width=1000,
            margin=dict(l=250, r=40, t=50, b=40),
            xaxis=dict(dtick=1, title="Number of Enriched Genes"),
            yaxis_title=None,
        )
        st.plotly_chart(fig, use_container_width=False)

    # --- Cross-Module Comparison ---
    st.header("Cross-Module Enrichment Comparison")
    database_to_compare = st.selectbox("Select Pathway Database", enrich_df['Source'].unique())
    enrich_compare_df = enrich_df[enrich_df['Source'] == database_to_compare]

    grouped_cmp = (
        enrich_compare_df
        .groupby(['Module', 'Path_Description'])['Gene']
        .nunique()
        .reset_index(name="Enriched_Genes")
    )

    pivot_cmp = grouped_cmp.pivot(index="Path_Description", columns="Module", values="Enriched_Genes").fillna(0)
    max_pathways = enrich_compare_df['Path_Description'].nunique()
    if max_pathways > 1:
        top_n = st.slider("Number of top pathways to show", 1, max_pathways, max(1, max_pathways))
    else:
        top_n=1

    pivot_cmp['Total'] = pivot_cmp.sum(axis=1)
    pivot_cmp = pivot_cmp.sort_values("Total", ascending=False).drop(columns="Total").head(top_n)




    # --- Grouped Bar Chart ---
    st.subheader("Grouped Bar Chart")
    fig_cmp = go.Figure()

    # Use Crest palette for modules
    modules = pivot_cmp.columns.tolist()
    colors = px.colors.sequential.Blues[1:] # GnBu alternatively
    module_colors = {module: colors[i % len(colors)] for i, module in enumerate(modules)}

    for module_name in modules:
        fig_cmp.add_trace(go.Bar(
            x=pivot_cmp.index,
            y=pivot_cmp[module_name],
            name=str(module_name),
            marker_color=module_colors[module_name],  # set color from palette
            hovertemplate=f"Module: {module_name}<br>Pathway: %{{x}}<br>Gene Count: %{{y}}<extra></extra>"
        ))

    fig_cmp.update_layout(
        legend_title_text='Module',
        barmode='group',
        title=f"Top {top_n} Enriched Pathways Across Modules - {database_to_compare}",
        xaxis_title="Pathway",
        yaxis_title="Number of Enriched Genes",
        xaxis_tickangle=-45,
        height=10 * top_n + 600,
        margin=dict(l=50, r=40, t=50, b=200),
    )
    st.plotly_chart(fig_cmp, use_container_width=True)


    # --- Stacked Bar Chart ---
    st.subheader("Stacked Bar Chart")
    fig_stacked = go.Figure()
    for module_name in modules:
        fig_stacked.add_trace(go.Bar(
            x=pivot_cmp.index,
            y=pivot_cmp[module_name],
            name=str(module_name),
            marker_color=module_colors[module_name],
            hovertemplate=f"Module: {module_name}<br>Pathway: %{{x}}<br>Gene Count: %{{y}}<extra></extra>"
        ))

    fig_stacked.update_layout(
        legend_title_text='Module',
        barmode='stack',
        title=f"Top {top_n} Enriched Pathways - {database_to_compare}",
        xaxis_title="Pathway",
        yaxis_title="Number of Enriched Genes",
        xaxis_tickangle=-45,
        height=10 * top_n + 600,
        margin=dict(l=50, r=40, t=50, b=200),
    )
    st.plotly_chart(fig_stacked, use_container_width=True)

    # --- Dot Plot: Modules vs Pathways (dot size proportional to enriched genes) ---
    st.subheader("Dot Plot: Number of Genes per Module-Pathway")
    # Use a color palette for pathways
    pathways = pivot_cmp.index.tolist()
    pathway_colors = px.colors.sequential.haline[:len(pathways)]  # Use a sequential palette
    color_map = {pathway: pathway_colors[i % len(pathway_colors)] for i, pathway in enumerate(pathways)}

    # Determine dot size scaling
    max_count = pivot_cmp.max().max()
    min_marker_size = 10
    max_marker_size = 40

    fig_dot = go.Figure()
    pathway_order = pivot_cmp.sum(axis=1).sort_values(ascending=False).index.tolist()

    for module_name in pivot_cmp.columns:
        for pathway, count in pivot_cmp[module_name].items():
            if count > 0:  # <-- only add dots for non-zero counts
                # Scale marker size proportional to count
                size = min_marker_size + (max_marker_size - min_marker_size) * (count / max_count)
                fig_dot.add_trace(go.Scatter(
                    x=[module_name],
                    y=[pathway],
                    mode='markers',
                    marker=dict(size=size, color=color_map[pathway], opacity=0.7),
                    hovertemplate=f"Module: {module_name}<br>Pathway: {pathway}<br>Enriched Genes: {int(count)}<extra></extra>",
                    showlegend=False
                ))

        # seed the y-axis categories so order is respected
    fig_dot.add_trace(go.Scatter(
        x=[pivot_cmp.columns[0]] * len(pathway_order),
        y=pathway_order,
        mode="markers",
        marker=dict(size=[0]*len(pathway_order), opacity=0),
        hoverinfo="skip",
        showlegend=False
))
    # Lock y-axis order: most → least at the top
    fig_dot.update_yaxes(
        categoryorder="array",
        categoryarray=pathway_order,
        autorange="reversed"  # puts the first category at the top
    )

    fig_dot.update_layout(
        xaxis_title="Module",
        yaxis_title="Pathway",
        yaxis=dict(autorange="reversed"),  # top pathways on top
        title=f"Dot Plot of Genes per Module-Pathway ({database_to_compare})",
        height=10 * top_n + 600,
        margin=dict(l=200, r=40, t=50, b=50),
    )

    st.plotly_chart(fig_dot, use_container_width=True)
    

    # --- Heatmap ---
    palette = px.colors.sequential.haline[::-1]  # reverse the haline palette
    custom = ["white"] + palette

    st.subheader("Heatmap: Modules vs Pathways")
    fig_heatmap = go.Figure(data=go.Heatmap(
        z=pivot_cmp.values,      # keep the original orientation
        x=pivot_cmp.columns,     # modules on x-axis
        y=pivot_cmp.index,       # pathways on y-axis
        colorscale=custom,
        colorbar=dict(title="Enriched Genes"),
        hovertemplate="Module: %{x}<br>Pathway: %{y}<br>Genes: %{z}<extra></extra>"
    ))
    fig_heatmap.update_layout(
        title=f"Heatmap of Enriched Genes - Top {top_n} Pathways ({database_to_compare})",
        xaxis_title="Module",
        yaxis_title="Pathway",
        yaxis=dict(autorange="reversed"),
        height=700 + 20 * len(pivot_cmp.index),
        margin=dict(l=200, r=40, t=50, b=200),  # increase left margin for long pathway names
    )
    st.plotly_chart(fig_heatmap, use_container_width=True)


