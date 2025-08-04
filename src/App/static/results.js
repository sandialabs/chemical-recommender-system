// © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
// SPDX-License-Identifier: BSD-3-Clause

document.addEventListener("DOMContentLoaded", function() {
  // Get embedded data from the template
  const finalresultsElement = document.getElementById("finalresults-data");
  const metadataElement = document.getElementById("results-metadata");
  const queryElement = document.getElementById("query-data");
  
  if (!finalresultsElement || !metadataElement || !queryElement) {
    console.error("Required data elements not found in template");
    return;
  }

  try {
    console.log("[DEBUG] === FRONTEND DATA RECEIVED ===");
    console.log("[DEBUG] finalresultsElement textContent:", finalresultsElement.textContent);
    console.log("[DEBUG] metadataElement textContent:", metadataElement.textContent);
    console.log("[DEBUG] queryElement textContent:", queryElement.textContent);
    
    const finalResults = JSON.parse(finalresultsElement.textContent);
    const metadata = JSON.parse(metadataElement.textContent);
    const queryData = JSON.parse(queryElement.textContent);
    
    console.log("[DEBUG] Final results parsed:", finalResults);
    console.log("[DEBUG] Final results type:", typeof finalResults);
    console.log("[DEBUG] Final results length:", finalResults ? finalResults.length : 'null/undefined');
    console.log("[DEBUG] Metadata:", metadata);
    console.log("[DEBUG] Metadata length:", metadata ? metadata.length : 'null/undefined');
    console.log("[DEBUG] Query data:", queryData);

    if (finalResults && finalResults.length > 0) {
        console.log("[DEBUG] First result:", finalResults[0]);
        console.log("[DEBUG] First result type:", typeof finalResults[0]);
        if (typeof finalResults[0] === 'object') {
            console.log("[DEBUG] First result keys:", Object.keys(finalResults[0]));
        }
    } else {
        console.log("[DEBUG] No results in finalResults array");
    }
    
    // Create cards for the chemical results
    createCards(finalResults, metadata, queryData);
    
    // Create query card
    createQCard(queryData);
    
  } catch (error) {
    console.error("Error parsing embedded data:", error);
  }
});

function createQCard(queryData) {
  const cardContainer = document.getElementById("card-container");
  if (!cardContainer) {
    console.error("card-container element not found");
    return;
  }
  
  // Clear existing cards (but keep the header)
  const existingCards = cardContainer.querySelectorAll('.card');
  existingCards.forEach(card => card.remove());
  
  const cardDiv = document.createElement("div");
  cardDiv.classList.add("card");
  
  console.log("[DEBUG] Query data received:", queryData);
  console.log("[DEBUG] qcid:", queryData.qcid);
  console.log("[DEBUG] qcid type:", typeof queryData.qcid);
  
  if (queryData.qcid && queryData.qcid !== -1 && queryData.qcid !== "-1") {
    // Clean the CID to ensure it's an integer
    let cleanCid = queryData.qcid;
    try {
      if (typeof cleanCid === 'string' && cleanCid.includes('.')) {
        cleanCid = parseInt(parseFloat(cleanCid));
      } else if (typeof cleanCid === 'number') {
        cleanCid = parseInt(cleanCid);
      }
    } catch (e) {
      console.log("[DEBUG] CID conversion error:", e);
    }
    
    const name = document.createElement("h3");
    name.innerHTML = "<strong>Query Compound</strong>";
    
    const cid = document.createElement("p");
    cid.innerHTML = '<strong>CID: </strong> <a href="https://pubchem.ncbi.nlm.nih.gov/compound/' + 
                    cleanCid + '" class="compound-id">' + cleanCid + "</a>";
    
    const image = document.createElement("img");
    image.classList.add("image");
    image.src = "https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid=" + cleanCid + "&t=l";
    image.alt = "Query Chemical Structure";
    
    cardDiv.appendChild(name);
    cardDiv.appendChild(cid);
    cardDiv.appendChild(image);
  } else {
    const name = document.createElement("h3");
    name.innerHTML = "<strong>Query: </strong>" + (queryData.queryval || "Unknown");
    console.log("[DEBUG] Query compound not found in PubChem, using query value:", queryData.queryval);
    const cid = document.createElement("p");
    cid.innerHTML = "<strong>CID: </strong> Not Found in PubChem";
    fetch('/api/query_image')
      .then(response => response.json())
      .then(data => {
        if (data && data.image) {
          console.log("[DEBUG] Query image data URL received:", data.image);
          const img = document.createElement("img");
          img.classList.add("image");
          img.src = data.image; // Use the data URL directly
          img.alt = "Query Chemical Structure";
          cardDiv.appendChild(img);
        }
      })
      .catch(error => {
        console.error("Error fetching query image bytes:", error);
      });
    cardDiv.appendChild(name);
    cardDiv.appendChild(cid);
  }
  
  cardContainer.appendChild(cardDiv);
}

function createCards(finalResults, metadata, queryData) {
  const cardsContainer = document.getElementById("cards-container");
  if (!cardsContainer) {
    console.error("cards-container element not found");
    return;
  }
  
  // Clear existing cards (but keep the header)
  const existingCards = cardsContainer.querySelectorAll('.card');
  existingCards.forEach(card => card.remove());
  
  if (!finalResults || !Array.isArray(finalResults) || finalResults.length === 0) {
    console.log("No results to display in cards");
    return;
  }
  
  // Get the query CID for comparison
  const queryCid = queryData && queryData.qcid ? String(queryData.qcid) : null;
  console.log("[DEBUG] Query CID for filtering:", queryCid);
  
  // Filter out the query compound from the results
  const results = finalResults.filter((item, index) => {
    // Skip if it's the query CID or has "Query" as the CID
    if (item.cid) {
      // Handle "Query" string
      if (item.cid === "Query") {
        console.log("[DEBUG] Filtering out 'Query' compound:", item.cid);
        return false;
      }
      
      // Handle CID comparison (normalize both to strings for comparison)
      if (queryCid) {
        const itemCidStr = String(item.cid).replace('.0', ''); // Remove .0 if present
        const queryCidStr = String(queryCid).replace('.0', '');
        if (itemCidStr === queryCidStr) {
          console.log("[DEBUG] Filtering out query compound by CID:", item.cid, "matches", queryCid);
          return false;
        }
      }
    }
    return true;
  });
  
  console.log("[DEBUG] Filtered results for cards:", results.length, "from original:", finalResults.length);
  
  results.forEach((item, index) => {
    const cardDiv = document.createElement("div");
    cardDiv.classList.add("card");

    // Add rank number
    const rank = document.createElement("h2");
    rank.innerHTML = '<strong style="color:#043885;">' + 
                     (index + 1) + ". </strong>";
    rank.style.marginBottom = "10px";
    rank.style.marginTop = "0px";
    cardDiv.appendChild(rank);

    // Extract CID from the result (now plain CID, no HTML formatting)
    let cidValue = item.cid;
    if (typeof cidValue === 'number') {
      cidValue = String(parseInt(cidValue));
    } else if (typeof cidValue === 'string') {
      try {
        if (cidValue.includes('.') && !isNaN(parseFloat(cidValue))) {
          cidValue = String(parseInt(parseFloat(cidValue)));
        }
      } catch (e) {}
    }

    // Add CID with link to PubChem
    const cid = document.createElement("p");
    cid.innerHTML = '<strong>PubChem CID: </strong> <a href="https://pubchem.ncbi.nlm.nih.gov/compound/' + 
                    cidValue + '" class="compound-id">' + cidValue + "</a>";
    cardDiv.appendChild(cid);

    // Add chemical structure image
    const image = document.createElement("img");
    image.classList.add("image");
    image.src = "https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid=" + cidValue + "&t=l";
    image.alt = "Chemical Structure";
    cardDiv.appendChild(image);

    // Add a divider
    const divider = document.createElement("div");
    divider.classList.add("divider");
    cardDiv.appendChild(divider);

    // Add other properties (but limit to most important ones for cards)
    const importantFields = ['final_score', 'structural_similarity', 'mw_similarity', 'thermo_similarity'];
    metadata.forEach((field) => {
      if (importantFields.includes(field.key) && item[field.key] !== undefined && item[field.key] !== null && item[field.key] !== "") {
        const elem = document.createElement("p");
        let value = item[field.key];
        if (field.type === "number" && typeof value === "number") {
          value = value.toFixed(3);
        }
        elem.innerHTML = '<strong>' + field.label + ': </strong><span class="score-value">' + value + '</span>';
        elem.classList.add("similarity-score");
        cardDiv.appendChild(elem);
      }
    });

    // Make card clickable for modal popup
    cardDiv.style.cursor = "pointer";
    cardDiv.addEventListener("click", function() {
      showCompoundModal(item, metadata);
    });

    cardsContainer.appendChild(cardDiv);
  });
}

// Modal logic
function showCompoundModal(item, metadata) {
  const modal = document.getElementById("compound-modal");
  const modalBody = document.getElementById("modal-body");
  if (!modal || !modalBody) return;

  // Clear previous content
  modalBody.innerHTML = "";

  // Create two columns: left for table, right for images
  const detailsCol = document.createElement("div");
  detailsCol.className = "modal-details-col";
  const imgCol = document.createElement("div");
  imgCol.className = "modal-img-col";


  // Structure image
  let cidValue = item.cid;
  if (typeof cidValue === 'number') cidValue = String(parseInt(cidValue));
  else if (typeof cidValue === 'string') {
    try {
      if (cidValue.includes('.') && !isNaN(parseFloat(cidValue))) {
        cidValue = String(parseInt(parseFloat(cidValue)));
      }
    } catch (e) {}
  }
  const img = document.createElement("img");
  img.src = "https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid=" + cidValue + "&t=l";
  img.alt = "Chemical Structure";
  img.className = "modal-structure-img";
  imgCol.appendChild(img);

  // Similarity map image (fetched from backend)
  const simMapDiv = document.createElement("div");
  simMapDiv.style.width = "100%";
  simMapDiv.style.display = "flex";
  simMapDiv.style.justifyContent = "center";
  simMapDiv.style.alignItems = "center";
  simMapDiv.style.minHeight = "110px";
  simMapDiv.style.background = "#f6f8fa";
  simMapDiv.style.borderRadius = "8px";
  simMapDiv.style.marginTop = "8px";
  simMapDiv.innerHTML = '<span style="color:#b0b8c9;font-size:1.05em;font-style:italic;">Loading similarity map...</span>';
  imgCol.appendChild(simMapDiv);

  // Fetch similarity map image from backend
  fetch(`/api/similarity_map?cid=${encodeURIComponent(cidValue)}`)
    .then(resp => resp.json())
    .then(data => {
      if (data && data.image) {
        simMapDiv.innerHTML = `<img src="${data.image}" alt="Similarity Map" style="height:100%;max-width:300px;border-radius:8px;box-shadow:0 2px 8px rgba(0,0,0,0.08);background:#fff;">`;
      } else if (data && data.error) {
        simMapDiv.innerHTML = `<span style="color:#d9534f;font-size:1.05em;font-style:italic;">Similarity map unavailable: ${data.error}</span>`;
      } else {
        simMapDiv.innerHTML = '<span style="color:#b0b8c9;font-size:1.05em;font-style:italic;">Similarity map unavailable</span>';
      }
    })
    .catch((err) => {
      simMapDiv.innerHTML = '<span style="color:#d9534f;font-size:1.05em;font-style:italic;">Similarity map unavailable (network error)</span>';
    });


  // Use a table for all available fields
  const table = document.createElement("table");
  table.className = "modal-details-table";
  metadata.forEach((field) => {
    if (item[field.key] !== undefined && item[field.key] !== null && item[field.key] !== "") {
      const tr = document.createElement("tr");
      const th = document.createElement("th");
      th.textContent = field.label;
      const td = document.createElement("td");
      td.textContent = (typeof item[field.key] === 'number' ? item[field.key].toFixed(4) : item[field.key]);
      tr.appendChild(th);
      tr.appendChild(td);
      table.appendChild(tr);
    }
  });
  detailsCol.appendChild(table);

  // Append detailsCol first (left), then imgCol (right)
  modalBody.appendChild(detailsCol);
  modalBody.appendChild(imgCol);

  // Show modal
  modal.style.display = "block";
}

// Close modal logic
document.addEventListener("DOMContentLoaded", function() {
  const modal = document.getElementById("compound-modal");
  const closeBtn = document.getElementById("close-modal");
  if (closeBtn && modal) {
    closeBtn.onclick = function() {
      modal.style.display = "none";
    };
    window.onclick = function(event) {
      if (event.target === modal) {
        modal.style.display = "none";
      }
    };
  }
});

