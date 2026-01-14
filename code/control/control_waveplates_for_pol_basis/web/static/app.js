
let currentPol = null;

function selectPol(pol) {
    currentPol = pol;

    document.querySelectorAll(".pol-btn").forEach(btn =>
        btn.classList.remove("active")
    );

    event.target.classList.add("active");

    document.getElementById("qwp").value = SETTINGS[pol].QWP;
    document.getElementById("hwp").value = SETTINGS[pol].HWP;

    fetch("/select", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ pol })
    });
}

function savePol() {
    if (!currentPol) return;

    SETTINGS[currentPol].QWP = document.getElementById("qwp").value;
    SETTINGS[currentPol].HWP = document.getElementById("hwp").value;

    fetch("/save", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
            pol: currentPol,
            QWP: SETTINGS[currentPol].QWP,
            HWP: SETTINGS[currentPol].HWP
        })
    }).then(() => flashButton());
}

function flashButton() {
    const btn = document.querySelector(".pol-btn.active");
    btn.classList.add("saved");
    setTimeout(() => btn.classList.remove("saved"), 300);
}

function exitApp() {
    fetch("/shutdown", { method: "POST" }).then(() => {
        window.close();
    });
}
function updateLogs() {
    fetch("/logs")
        .then(response => response.json())
        .then(messages => {
            const logDiv = document.getElementById("log");
            logDiv.textContent = messages.join("\n");
            logDiv.scrollTop = logDiv.scrollHeight; // auto-scroll
        })
        .catch(err => console.error(err));
}

// poll every 500 ms
setInterval(updateLogs, 500);
