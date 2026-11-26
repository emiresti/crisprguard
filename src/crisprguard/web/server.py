"""Zero-dependency HTTP server and JSON API router for the CrisprGuard Web Laboratory."""

from __future__ import annotations

from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
from typing import Any
from urllib.parse import urlparse
import webbrowser

from crisprguard import __version__
from crisprguard.design import (
    ENZYME_SPECS,
    design_guides,
    format_guides_csv,
    format_guides_fasta,
    predict_mmej,
)
from crisprguard.offtarget import scan_offtargets
from crisprguard.samples import list_samples

STATIC_DIR = Path(__file__).resolve().parent / "static"
MAX_BODY_BYTES = 2 * 1024 * 1024  # 2 MiB guard

MIME_TYPES: dict[str, str] = {
    ".html": "text/html; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".js": "application/javascript; charset=utf-8",
    ".svg": "image/svg+xml",
    ".json": "application/json; charset=utf-8",
    ".png": "image/png",
    ".ico": "image/x-icon",
}


def _json_response(status: int, payload: Any) -> tuple[int, bytes, str]:
    body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
    return status, body, "application/json; charset=utf-8"


def _parse_json_body(body: bytes) -> dict[str, Any]:
    if not body or not body.strip():
        raise ValueError("Request JSON body is required.")
    if len(body) > MAX_BODY_BYTES:
        raise ValueError("Request payload exceeds maximum allowed size (2 MiB).")
    try:
        parsed = json.loads(body.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"Malformed JSON payload: {exc}") from exc
    if not isinstance(parsed, dict):
        raise ValueError("JSON body must be an object.")
    return parsed


def _serve_static(rel_path: str) -> tuple[int, bytes, str]:
    clean_rel = "index.html" if rel_path in {"", "/"} else rel_path.lstrip("/")
    candidate = (STATIC_DIR / clean_rel).resolve()
    try:
        candidate.relative_to(STATIC_DIR)
    except ValueError:
        return _json_response(HTTPStatus.FORBIDDEN, {"error": "Forbidden path."})
    if not candidate.is_file():
        return _json_response(HTTPStatus.NOT_FOUND, {"error": f"Asset '{clean_rel}' not found."})
    mime = MIME_TYPES.get(candidate.suffix.lower(), "application/octet-stream")
    return HTTPStatus.OK, candidate.read_bytes(), mime


def route_request(
    method: str,
    raw_path: str,
    body: bytes = b"",
) -> tuple[int, bytes, str]:
    """Pure, testable HTTP router returning (status_code, response_bytes, content_type)."""
    verb = method.upper()
    path = urlparse(raw_path).path

    if verb == "OPTIONS":
        return HTTPStatus.NO_CONTENT, b"", "text/plain; charset=utf-8"

    if path == "/api/health":
        if verb != "GET":
            return _json_response(HTTPStatus.METHOD_NOT_ALLOWED, {"error": "Use GET."})
        return _json_response(
            HTTPStatus.OK,
            {
                "status": "ok",
                "name": "CrisprGuard",
                "version": __version__,
                "enzymes": {k: v.to_dict() for k, v in ENZYME_SPECS.items()},
            },
        )

    if path == "/api/samples":
        if verb != "GET":
            return _json_response(HTTPStatus.METHOD_NOT_ALLOWED, {"error": "Use GET."})
        return _json_response(
            HTTPStatus.OK,
            {"samples": [s.to_dict() for s in list_samples()]},
        )

    if path == "/api/design":
        if verb != "POST":
            return _json_response(HTTPStatus.METHOD_NOT_ALLOWED, {"error": "Use POST."})
        try:
            data = _parse_json_body(body)
            sequence = str(data.get("sequence", ""))
            if not sequence.strip():
                raise ValueError("Field 'sequence' is required and cannot be empty.")
            enzyme = str(data.get("enzyme", "spcas9"))
            background = data.get("background")
            bg_str = str(background) if background is not None else None
            max_mm = int(data.get("max_mismatches", 4))
            include_noncanonical = bool(data.get("include_noncanonical_pam", True))
            mmej_window = int(data.get("mmej_window", 25))
            min_gc = float(data.get("min_gc", 0.0))
            max_gc = float(data.get("max_gc", 100.0))
            top_n_raw = data.get("top_n")
            top_n = int(top_n_raw) if top_n_raw is not None else None
            target_name = data.get("target_name")

            report = design_guides(
                sequence=sequence,
                enzyme=enzyme,
                background=bg_str,
                max_mismatches=max_mm,
                include_noncanonical_pam=include_noncanonical,
                mmej_window=mmej_window,
                top_n=top_n,
                min_gc=min_gc,
                max_gc=max_gc,
                target_name=str(target_name) if target_name else None,
            )
            return _json_response(HTTPStatus.OK, report.to_dict())
        except (ValueError, TypeError) as exc:
            return _json_response(HTTPStatus.BAD_REQUEST, {"error": str(exc)})

    if path == "/api/offtarget":
        if verb != "POST":
            return _json_response(HTTPStatus.METHOD_NOT_ALLOWED, {"error": "Use POST."})
        try:
            data = _parse_json_body(body)
            guide = str(data.get("guide", ""))
            background = str(data.get("background", ""))
            if not guide.strip():
                raise ValueError("Field 'guide' is required.")
            if not background.strip():
                raise ValueError("Field 'background' is required.")
            enzyme = str(data.get("enzyme", "spcas9"))
            max_mm = int(data.get("max_mismatches", 4))
            include_noncanonical = bool(data.get("include_noncanonical_pam", True))
            top_n = int(data.get("top_n", 25))
            exclude_one_exact = bool(data.get("exclude_one_exact", False))

            summary = scan_offtargets(
                guide=guide,
                background=background,
                enzyme=enzyme,
                max_mismatches=max_mm,
                include_noncanonical_pam=include_noncanonical,
                exclude_one_exact=exclude_one_exact,
                top_n=top_n,
            )
            return _json_response(HTTPStatus.OK, summary.to_dict())
        except (ValueError, TypeError) as exc:
            return _json_response(HTTPStatus.BAD_REQUEST, {"error": str(exc)})

    if path == "/api/mmej":
        if verb != "POST":
            return _json_response(HTTPStatus.METHOD_NOT_ALLOWED, {"error": "Use POST."})
        try:
            data = _parse_json_body(body)
            sequence = str(data.get("sequence", ""))
            if not sequence.strip():
                raise ValueError("Field 'sequence' is required.")
            if "cut_pos" not in data or data["cut_pos"] is None:
                raise ValueError("Field 'cut_pos' is required.")
            cut_pos = int(data["cut_pos"])
            window = int(data.get("window", 25))
            top_n = int(data.get("top_n", 12))

            result = predict_mmej(sequence=sequence, cut_pos=cut_pos, window=window, top_n=top_n)
            return _json_response(HTTPStatus.OK, result.to_dict())
        except (ValueError, TypeError) as exc:
            return _json_response(HTTPStatus.BAD_REQUEST, {"error": str(exc)})

    if path == "/api/export":
        if verb != "POST":
            return _json_response(HTTPStatus.METHOD_NOT_ALLOWED, {"error": "Use POST."})
        try:
            data = _parse_json_body(body)
            sequence = str(data.get("sequence", ""))
            if not sequence.strip():
                raise ValueError("Field 'sequence' is required.")
            fmt = str(data.get("format", "csv")).strip().lower()
            if fmt not in {"csv", "fasta"}:
                raise ValueError("Field 'format' must be 'csv' or 'fasta'.")
            enzyme = str(data.get("enzyme", "spcas9"))
            background = data.get("background")
            bg_str = str(background) if background is not None else None
            max_mm = int(data.get("max_mismatches", 4))
            top_n_raw = data.get("top_n")
            top_n = int(top_n_raw) if top_n_raw is not None else None

            report = design_guides(
                sequence=sequence,
                enzyme=enzyme,
                background=bg_str,
                max_mismatches=max_mm,
                top_n=top_n,
            )
            if fmt == "csv":
                text = format_guides_csv(report.guides)
                return HTTPStatus.OK, text.encode("utf-8"), "text/csv; charset=utf-8"
            text = format_guides_fasta(report.guides, target_name=report.target_name)
            return HTTPStatus.OK, text.encode("utf-8"), "text/plain; charset=utf-8"
        except (ValueError, TypeError) as exc:
            return _json_response(HTTPStatus.BAD_REQUEST, {"error": str(exc)})

    if path.startswith("/api/"):
        return _json_response(HTTPStatus.NOT_FOUND, {"error": f"Unknown API endpoint '{path}'."})

    if verb == "GET":
        return _serve_static(path)

    return _json_response(HTTPStatus.METHOD_NOT_ALLOWED, {"error": "Method not allowed."})


class CrisprGuardRequestHandler(BaseHTTPRequestHandler):
    """HTTP handler delegating all requests to the pure `route_request` router."""

    def log_message(self, format: str, *args: Any) -> None:  # noqa: A002
        # Quiet clean server logging
        return

    def _dispatch(self, method: str) -> None:
        length = int(self.headers.get("Content-Length", "0") or "0")
        if length > MAX_BODY_BYTES:
            status, payload, ctype = _json_response(
                HTTPStatus.REQUEST_ENTITY_TOO_LARGE,
                {"error": "Request entity too large."},
            )
        else:
            body = self.rfile.read(length) if length > 0 else b""
            status, payload, ctype = route_request(method, self.path, body)

        self.send_response(int(status))
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        if payload:
            self.wfile.write(payload)

    def do_GET(self) -> None:  # noqa: N802
        self._dispatch("GET")

    def do_POST(self) -> None:  # noqa: N802
        self._dispatch("POST")

    def do_OPTIONS(self) -> None:  # noqa: N802
        self._dispatch("OPTIONS")


def serve(host: str = "127.0.0.1", port: int = 8765, open_browser: bool = True) -> None:
    """Start the local CrisprGuard interactive web laboratory server."""
    server = ThreadingHTTPServer((host, port), CrisprGuardRequestHandler)
    url = f"http://{host}:{port}"
    print(f"CrisprGuard Web Laboratory v{__version__} running at {url}")
    print("Press Ctrl+C to stop.")
    if open_browser:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down CrisprGuard server.")
    finally:
        server.server_close()
