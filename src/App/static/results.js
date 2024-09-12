// © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
// SPDX-License-Identifier: BSD-3-Clause

document.addEventListener("DOMContentLoaded", fetchAndCreateCards);
document.addEventListener("DOMContentLoaded", fetchAndCreateCard);

function fetchAndCreateCards() {
  const apiEndpoint = "/api/chemicals";
  fetch(apiEndpoint)
    .then((response) => response.json())
    .then((data) => {
      createCards(data);
      document.removeEventListener("DOMContentLoaded", fetchAndCreateCards);
    })
    .catch((error) => console.error("Error fetching data", error));
}

function fetchAndCreateCard() {
  const apiEndpoint = "/api/query";
  fetch(apiEndpoint)
    .then((response) => response.json())
    .then((data) => {
      console.log(data);
      if (data === null) {
        console.log("empty q card");
        createEmptyQCard();
      } else {
        createQCard(data);
      }
      document.removeEventListener("DOMContentLoaded", fetchAndCreateCard);
    })
    .catch((error) => console.error("Error fetching data", error));
}

function createEmptyQCard() {
  const cardContainer = document.getElementById("card-container");
  while (cardContainer.firstChild) {
    cardContainer.removeChild(cardContainer.firstChild);
  }
  cardContainer.innerHTML = "";

  const cardDiv = document.createElement("div");
  cardDiv.classList.add("card");

  const name = document.createElement("h3");
  name.innerHTML = "<strong>Name: </strong> Not Found in PubChem";

  const cid = document.createElement("p");
  cid.innerHTML = "<strong>CID: </strong> Not Found in PubChem";

  // cardDiv.appendChild(rank)
  cardDiv.appendChild(name);
  cardDiv.appendChild(cid);
  console.log(cardDiv);
  cardContainer.appendChild(cardDiv);
}

function createQCard(data) {
  const cardContainer = document.getElementById("card-container");
  while (cardContainer.firstChild) {
    cardContainer.removeChild(cardContainer.firstChild);
  }
  cardContainer.innerHTML = "";
  data.forEach((chemical) => {
    const cardDiv = document.createElement("div");
    cardDiv.classList.add("card");

    const name = document.createElement("h3");
    name.innerHTML = "<strong>Name: </strong>" + chemical.name;

    const cid = document.createElement("p");
    cid.innerHTML =
      '<strong>CID: </strong> <a href="https://pubchem.ncbi.nlm.nih.gov/compound/' +
      chemical.cid +
      '">' +
      chemical.cid +
      "</a>";

    const image = document.createElement("img");
    image.classList.add("image");
    image.src = chemical.im;
    image.alt = "Chemical Image";

    // cardDiv.appendChild(rank)
    cardDiv.appendChild(name);
    cardDiv.appendChild(cid);
    cardDiv.appendChild(image);

    cardContainer.appendChild(cardDiv);
  });
}

function createCards(data) {
  const cardsContainer = document.getElementById("cards-container");
  while (cardsContainer.firstChild) {
    cardsContainer.removeChild(cardsContainer.firstChild);
  }
  cardsContainer.innerHTML = "";
  var count = 0;
  data.forEach((chemical) => {
    count = count + 1;
    const cardDiv = document.createElement("div");
    cardDiv.classList.add("card");

    const rank = document.createElement("h2");
    rank.innerHTML =
      '<strong style="position:relative; right:11.9vw; color:#043885;">' +
      count +
      ". </strong>";
    rank.style.marginBottom = "0px";
    rank.style.marginTop = "5px";

    const name = document.createElement("h3");
    name.innerHTML = "<strong>Name: </strong>" + chemical.name;
    name.style.marginTop = "0px";

    const cid = document.createElement("p");
    cid.innerHTML =
      '<strong>CID: </strong> <a href="https://pubchem.ncbi.nlm.nih.gov/compound/' +
      chemical.cid +
      '">' +
      chemical.cid +
      "</a>";

    const image = document.createElement("img");
    image.classList.add("image");
    image.src = chemical.im;
    image.alt = "Chemical Image";

    cardDiv.appendChild(rank);
    cardDiv.appendChild(name);
    cardDiv.appendChild(cid);
    cardDiv.appendChild(image);

    cardsContainer.appendChild(cardDiv);
  });
}
